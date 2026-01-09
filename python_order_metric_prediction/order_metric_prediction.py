# Module to predict the outcome of the WWW Monte Carlo algorithm
# this code is taken and modified from the following link:
# https://docs.pytorch.org/tutorials/beginner/hyperparameter_tuning_tutorial.html

# Import packages
from functools import partial
import os
from pathlib import Path
import tempfile
import numpy as np
import random
import h5py
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset
from ray import tune
from ray import train
from ray.tune import Checkpoint, get_checkpoint
from ray.tune.schedulers import ASHAScheduler
import ray.cloudpickle as pickle
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.metrics import r2_score


class NeuralNetwork(nn.Module):
    """
    Neural network with variable number of layers and fixed number of neurons
    per hidden layer.
    """
    def __init__(self, nr_inputs=10, nr_outputs=3, nr_neurons=10, nr_layers=2):
        super(NeuralNetwork, self).__init__()
        self.fc_in = nn.Linear(nr_inputs, nr_neurons)
        self.fc_hidden = nn.Linear(nr_neurons, nr_neurons)
        self.fc_out = nn.Linear(nr_neurons, nr_outputs)
        self.nr_layers = nr_layers

    def forward(self, x):
        x = F.relu(self.fc_in(x))
        for _ in range(self.nr_layers - 1):
            x = F.relu(self.fc_hidden(x))
        x = self.fc_out(x)
        return x


class NeuralNetworkVariableNrNeurons(nn.Module):
    """
    Neural network with variable number of layers and variable number of
    neurons per hidden layer.
    """
    def __init__(self, nr_inputs=10, nr_outputs=3, nr_neurons_vec=[8, 6]):
        super(NeuralNetworkVariableNrNeurons, self).__init__()

        layers = [nn.Linear(nr_inputs, nr_neurons_vec[0]), nn.ReLU()]

        for i in range(len(nr_neurons_vec) - 1):
            layers.append(nn.Linear(nr_neurons_vec[i], nr_neurons_vec[i+1]))
            layers.append(nn.ReLU())

        layers.append(nn.Linear(nr_neurons_vec[-1], nr_outputs))

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)


class WeightedMSE(nn.Module):
    """Custom loss function that applies loss weights to the components of the 
    predicted and true vectors.
    """
    def __init__(self, weights):
        super(WeightedMSE, self).__init__()
        self.weights = weights

    def forward(self, input, target):
        # Compute the loss
        loss = torch.mean(self.weights * (input - target) ** 2)
        return loss


def loss_batch(model, loss_func, xb, yb, opt=None):
    """
    Compute the loss for a batch of data. If an optimizer is provided, perform
    a backward pass and update the model parameters.
    """
    loss = loss_func(model(xb), yb)

    if opt is not None:
        loss.backward()
        opt.step()
        opt.zero_grad()

    return loss.item(), len(xb)


def fit(epochs, model, loss_func, opt, train_dl, valid_dl):
    """ Fit the model to the training data and evaluate on the validation data.
    """
    for epoch in range(epochs):
        model.train()
        for xb, yb in train_dl:
            loss_batch(model, loss_func, xb, yb, opt)

        model.eval()
        with torch.no_grad():
            losses, nums = zip(
                *[loss_batch(model, loss_func, xb, yb) for xb, yb in valid_dl]
            )
        val_loss = np.sum(np.multiply(losses, nums)) / np.sum(nums)

        print(epoch, val_loss)


def create_dataset_metric_prediction(
        load_data_path="../data/analysis_data/ctn/",
        load_filename="all_order_metrics.h5",
        neural_network_data_path="../data/neural_networks/ctn/",
        save_filename="order_metric_dataset.h5",
        test_fraction = 0.15,
        apply_pca = True,
        nr_pca_components=10,
        seed=42,
        order_metrics_vec = [
        "bond_length_std_vec",
        "bond_angle_std_vec",
        "dihedral_angle_entropy_vec",
        "bond_orientation_entropy_vec",
        "coordination_nr_mean_vec",
        "coordination_nr_std_vec",
        "q_l_mat_values",
        "q_l_mat_uncertainties",
        "vertex_homogeneity_metric_vec",
        "uncoordinated_neighbor_distance_vec",
        "ring_size_mean_vec",
        "ring_size_std_vec",
        "ring_radius_mean_vec",
        "ring_radius_std_vec",
        "critical_pore_radius_vec",
        "anisotropy_metric_from_structure_factor_vec",
        "anisotropy_metric_from_structure_factor_bonds_vec",
        "hyperuniformity_alpha_vec"]
        ):
    """ Create a dataset from an order metric dictionary as it is created by
    the Julia code.
    """

    # Load the data from the HDF5 file
    with h5py.File(load_data_path+load_filename, 'r') as hf:
        bond_bending_const_vec = np.array(hf["bond_bending_const_vec"])
        t_max_vec = np.array(hf["t_max_vec"])
        t_gradient_vec = np.array(hf["t_gradient_vec"])

        x_data = np.column_stack((bond_bending_const_vec, t_max_vec,
                                  t_gradient_vec))
        
        count = 0

        for order_metric in order_metrics_vec:

            # Order metrics related to the Steinhardt parameters are matrices
            # and not vectors and are therefore treated separately
            if (order_metric == "q_l_mat_values" 
                or order_metric == "q_l_mat_uncertainties"):
                current_y_data = np.array(hf[order_metric])
                nr_q_l = current_y_data.shape[1]

                current_order_metrics_names = [
                    f"{order_metric}_{i}" for i in range(nr_q_l)]
                
                if count == 0:
                    y_data = current_y_data
                    order_metrics_names = current_order_metrics_names
                else:
                    y_data = np.concatenate((y_data, current_y_data), axis=1)
                    order_metrics_names.extend(current_order_metrics_names)
                
                count += nr_q_l

            else:
                current_y_data = np.array(hf[order_metric]).reshape(-1, 1)

                if count == 0:
                    y_data = current_y_data
                    order_metrics_names = [order_metric]
                else:
                    y_data = np.concatenate((y_data, current_y_data), axis=1)
                    order_metrics_names.append(order_metric)
                count += 1
            
    # convert the data to float32
    x_data = x_data.astype(np.float32)
    y_data = y_data.astype(np.float32)

    nr_samples = y_data.shape[0]

    # shuffle the rows of the x_data and y_data arrays in the same way
    # to ensure that the data is still aligned
    rng = np.random.default_rng(seed)
    perm = rng.permutation(nr_samples)
    x_data = x_data[perm]
    y_data = y_data[perm]

    print(order_metrics_names)
    print(y_data.shape)
    print(x_data.shape)

    # Split the data into training and validation and test sets
    test_size = int(test_fraction * nr_samples)
    rest_size = nr_samples - test_size
    train_size = int(0.8 * rest_size)
    valid_size = nr_samples - train_size - test_size

    y_train, y_valid, y_test = (y_data[:train_size], 
                                y_data[train_size:train_size+valid_size], 
                                y_data[train_size+valid_size:])
    x_train, x_valid, x_test = (x_data[:train_size], 
                                x_data[train_size:train_size+valid_size], 
                                x_data[train_size+valid_size:])

    if apply_pca:
        # Apply scaler and PCA to the y_data
        y_train, y_valid, y_test = apply_pca_to_data(
            y_train, y_valid, y_test,
            save_filename, 
            neural_network_data_path=neural_network_data_path,
            nr_pca_components=nr_pca_components)
        
        # Modify save_filename to include PCA information
        save_filename = f"{save_filename[:-3]}_pca_{nr_pca_components}.h5"
    
    # save the data to an hdf5 file
    with h5py.File(neural_network_data_path+save_filename, 'w') as hf:
        # Create a dataset in the file and write the data to it
        dset = hf.create_dataset("y_train", data=y_train)
        dset = hf.create_dataset("y_valid", data=y_valid)
        dset = hf.create_dataset("y_test", data=y_test)
        dset = hf.create_dataset("x_train", data=x_train)
        dset = hf.create_dataset("x_valid", data=x_valid)
        dset = hf.create_dataset("x_test", data=x_test)
        dset = hf.create_dataset("order_metrics_names", 
            data=order_metrics_names)
        
    return


def create_dataset_single_file(
        load_data_path="../data/analysis_data/ctn/",
        load_filename="sample_name_order_metrics.h5",
        order_metrics = [
        "bond_length_std",
        "bond_angle_std",
        "dihedral_angle_entropy",
        "bond_orientation_entropy",
        "coordination_nr_mean",
        "coordination_nr_std",
        "q_l_vec_values",
        "q_l_vec_uncertainties",
        "vertex_homogeneity_metric",
        "uncoordinated_neighbor_distance",
        "ring_size_mean",
        "ring_size_std",
        "ring_radius_mean",
        "ring_radius_std",
        "critical_pore_radius",
        "anisotropy_metric_from_structure_factor",
        "anisotropy_metric_from_structure_factor_bonds",
        "hyperuniformity_alpha"]
        ):
    """ Create a dataset from an order metric dictionary for a single network 
    as it is created by the Julia code.
    """

    # Load the data from the HDF5 file
    with h5py.File(load_data_path+load_filename, 'r') as hf:
        
        count = 0

        for order_metric in order_metrics:

            # Order metrics related to the Steinhardt parameters are vectors
            # and not scalars and are therefore treated separately
            if (order_metric == "q_l_vec_values" 
                or order_metric == "q_l_vec_uncertainties"):
                current_y_data = np.array(hf["q_l_vec_values"])
                nr_q_l = current_y_data.shape[0]

                current_order_metrics_names = [
                    f"{order_metric}_{i}" for i in range(nr_q_l)]
                
                if count == 0:
                    y_data = current_y_data
                    order_metrics_names = current_order_metrics_names
                else:
                    y_data = np.concatenate((y_data, current_y_data), axis=0)
                    order_metrics_names.extend(current_order_metrics_names)
                
                count += nr_q_l

            else:
                current_y_data = np.array([float(hf[order_metric][()])])

                if count == 0:
                    y_data = current_y_data
                    order_metrics_names = [order_metric]
                else:
                    y_data = np.concatenate((y_data, current_y_data), axis=0)
                    order_metrics_names.append(order_metric)
                count += 1
    
    # convert the data to float32
    y_data = y_data.astype(np.float32)

    print(order_metrics_names)
    print(y_data.shape)

    return y_data


def apply_pca_to_data(y_train, y_valid, y_test, 
                      filename, 
                      neural_network_data_path="../data/neural_networks/ctn/",
                      nr_pca_components=10):
    """ Standardize y data and then apply PCA to reduce its dimension. Then 
    save the the scaler and the PCA transformation to files.
    """

    # Standardize output before PCA
    scaler = StandardScaler()
    y_train_scaled = scaler.fit_transform(y_train)
    y_valid_scaled = scaler.transform(y_valid)
    y_test_scaled = scaler.transform(y_test)

    # Apply PCA
    pca = PCA(n_components=nr_pca_components)
    y_train_pca = pca.fit_transform(y_train_scaled)
    y_valid_pca = pca.transform(y_valid_scaled)
    y_test_pca = pca.transform(y_test_scaled)

    print(f"Explained variance ratio: {pca.explained_variance_ratio_}")
    print(f"Sum of explained variance ratio: {np.sum(pca.explained_variance_ratio_)}")

    # Save both scaler and PCA
    scaler_filename = f"{filename[:-3]}_pca_{nr_pca_components}_scaler.pkl"
    with open(neural_network_data_path+scaler_filename, 'wb') as f:
        pickle.dump(scaler, f)

    pca_filename = f"{filename[:-3]}_pca_{nr_pca_components}.pkl"
    with open(neural_network_data_path+pca_filename, 'wb') as f:
        pickle.dump(pca, f)

    save_pca_to_h5(pca, scaler, filename, 
                   neural_network_data_path=neural_network_data_path,
                   nr_pca_components=nr_pca_components)

    return y_train_pca, y_valid_pca, y_test_pca


def save_pca_to_h5(
        pca,
        scaler,
        filename,
        neural_network_data_path="../data/neural_networks/ctn/",
        nr_pca_components=10):
    """ Save the PCA transformation and scaler to an HDF5 file."""

    pca_filename = f"{filename[:-3]}_pca_{nr_pca_components}_information.h5"
    with h5py.File(neural_network_data_path+pca_filename, 'w') as hf:
        # Save PCA components
        hf.create_dataset("components_", data=pca.components_)
        hf.create_dataset("explained_variance_", data=pca.explained_variance_)
        hf.create_dataset("explained_variance_ratio_", 
                          data=pca.explained_variance_ratio_)
        hf.create_dataset("mean_", data=pca.mean_)
        # Save scaler parameters
        hf.create_dataset("scale_", data=scaler.scale_)
        hf.create_dataset("mean_scaler_", data=scaler.mean_)
        
    return


def load_scaler_and_pca(
        filename,
        neural_network_data_path="../data/neural_networks/ctn/"):
    """ Load the scaler and PCA transformation."""

    scaler_filename = f"{filename[:-3]}_scaler.pkl"
    with open(neural_network_data_path+scaler_filename, 'rb') as f:
        scaler = pickle.load(f)
    pca_filename = f"{filename[:-3]}.pkl"
    with open(neural_network_data_path+pca_filename, 'rb') as f:
        pca = pickle.load(f)

    return scaler, pca


def load_dataset(
        neural_network_data_path="../data/neural_networks/ctn/", 
        filename="order_metric_dataset.h5"):
    """ Load the dataset for order metric prediction from an HDF5 file and
    convert it to PyTorch tensors and TensorDatasets.
    ️"""

    # check if file exists
    if Path(neural_network_data_path+filename).exists():
        print("File found.")
        full_path = neural_network_data_path + filename
    else:
        SCRIPT_DIR = Path(__file__).resolve().parent
        rel_or_abs=Path(neural_network_data_path)
        full_path = (rel_or_abs if rel_or_abs.is_absolute() 
                     else SCRIPT_DIR / rel_or_abs) / filename
        full_path = full_path.resolve(strict=False)
    
    # Load the data from the HDF5 file
    with h5py.File(full_path, 'r') as hf:
        x_train = hf["x_train"][:]
        x_valid = hf["x_valid"][:]
        x_test = hf["x_test"][:]
        y_train = hf["y_train"][:]
        y_valid = hf["y_valid"][:]
        y_test = hf["y_test"][:]

    # Convert the data to PyTorch tensors
    x_train, y_train, x_valid, y_valid, x_test, y_test = map(
        torch.tensor, (x_train, y_train, x_valid, y_valid, x_test, y_test))
    train_ds = TensorDataset(x_train, y_train)
    valid_ds = TensorDataset(x_valid, y_valid)
    test_ds = TensorDataset(x_test, y_test)

    return train_ds, valid_ds, test_ds


def get_weights_balanced():
    """ Set balanced weights for the loss function. The weights are normalized
    to sum to 1 and are chosen such that different kinds of order metrics are
    weighted equally:
    - bond_length_std_vec 15/100
    - bond_angle_std_vec 15/100
    - dihedral_angle_entropy_vec 5/100
    - q_l_mat_values 5/100
    - q_l_mat_uncertainties 5/100
    Metrics describing homogeneity: 
    - critical_pore_radius_vec 10/100
    - vertex_homogeneity_metric_vec 5/100
    - uncoordinated_neighbor_distance_vec 5/100
    - hyperuniformity_alpha_vec 2/100
    Metrics describing isotropy: 
    - bond_orientation_entropy_vec 10/100
    - anisotropy_metric_from_structure_factor_vec 5/100
    - anisotropy_metric_from_structure_factor_bonds_vec 5/100
    Metrics describing topology: 
    - ring_size_mean_vec 0.5/100
    - ring_size_std_vec 0.5/100
    - ring_radius_mean_vec 5/100
    - ring_radius_std_vec 5/100
    - coordination_nr_mean_vec 1/100
    - coordination_nr_std_vec 1/100
    """

    order_metrics_vec = [
        "bond_length_std_vec", # 0.15
        "bond_angle_std_vec", # 0.15
        "dihedral_angle_entropy_vec", # 0.05
        "bond_orientation_entropy_vec", # 0.10
        "coordination_nr_mean_vec", # 0.01
        "coordination_nr_std_vec", # 0.01
        "q_l_mat_values", # 0.05 * 1/13
        "q_l_mat_uncertainties", # 0.05 * 1/13
        "vertex_homogeneity_metric_vec", # 0.05
        "uncoordinated_neighbor_distance_vec", # 0.05
        "ring_size_mean_vec", # 0.005
        "ring_size_std_vec", # 0.005
        "ring_radius_mean_vec", # 0.05
        "ring_radius_std_vec", # 0.05
        "critical_pore_radius_vec", # 0.10
        "anisotropy_metric_from_structure_factor_vec", # 0.05
        "anisotropy_metric_from_structure_factor_bonds_vec", # 0.05
        "hyperuniformity_alpha_vec", # 0.02
        ]

    weights = torch.tensor([
        0.15, # bond_length_std_vec
        0.15, # bond_angle_std_vec
        0.05, # dihedral_angle_entropy_vec
        0.10, # bond_orientation_entropy_vec
        0.01, # coordination_nr_mean_vec
        0.01, # coordination_nr_std_vec
        ] + [0.05 * 1/13 for _ in range(13) # q_l_mat_values
        ] + [0.05 * 1/13 for _ in range(13) # q_l_mat_uncertainties
        ] + [
        0.05, # vertex_homogeneity_metric_vec
        0.05, # uncoordinated_neighbor_distance_vec
        0.005, # ring_size_mean_vec
        0.005, # ring_size_std_vec
        0.05, # ring_radius_mean_vec
        0.05, # ring_radius_std_vec
        0.10, # critical_pore_radius_vec
        0.05, # anisotropy_metric_from_structure_factor_vec
        0.05, # anisotropy_metric_from_structure_factor_bonds_vec
        0.02, # hyperuniformity_alpha_vec
        ])
    
    return weights


def get_weights_focussed():
    """ Set weights for the loss function that focus on bond lengths and 
    angles. The weights are normalized to sum to 1.
    Metrics describing the primitives: 
    - bond_length_std_vec 25/100
    - bond_angle_std_vec 25/100
    - dihedral_angle_entropy_vec 2/100
    - q_l_mat_values 2/100
    - q_l_mat_uncertainties 2/100
    Metrics describing homogeneity: 
    - vertex_homogeneity_metric_vec 6/100
    - uncoordinated_neighbor_distance_vec 6/100
    - critical_pore_radius_vec 7/100
    Metrics describing anisotropy: 
    - bond_orientation_entropy_vec 9/100
    - anisotropy_metric_from_structure_factor_vec 4.5/100
    - anisotropy_metric_from_structure_factor_bonds_vec 4.5/100
    Metrics describing hyperuniformity: 
    - hyperuniformity_alpha_vec 2/100
    Metrics describing ring statistics: 
    - ring_size_mean_vec 0.5/100
    - ring_size_std_vec 0.5/100
    - ring_radius_mean_vec 1/100
    - ring_radius_std_vec 1/100
    Coordination number: 
    - coordination_nr_mean_vec 1/100
    - coordination_nr_std_vec 1/100
    """

    order_metrics_vec = [
        "bond_length_std_vec", # 25/100
        "bond_angle_std_vec", # 25/100
        "dihedral_angle_entropy_vec", # 2/100
        "bond_orientation_entropy_vec", # 9/100
        "coordination_nr_mean_vec", # 1/100
        "coordination_nr_std_vec", # 1/100
        "q_l_mat_values", # 2/100 * 1/13
        "q_l_mat_uncertainties", # 2/100 * 1/13
        "vertex_homogeneity_metric_vec", # 6/100
        "uncoordinated_neighbor_distance_vec", # 6/100
        "ring_size_mean_vec", # 0.5/100
        "ring_size_std_vec", # 0.5/100
        "ring_radius_mean_vec", # 1/100
        "ring_radius_std_vec", # 1/100
        "critical_pore_radius_vec", # 7/100
        "anisotropy_metric_from_structure_factor_vec", # 4.5/100
        "anisotropy_metric_from_structure_factor_bonds_vec", # 4.5/100
        "hyperuniformity_alpha_vec", # 2/100
        ]

    weights = torch.tensor([
        25/100, # bond_length_std_vec
        25/100, # bond_angle_std_vec
        2/100, # dihedral_angle_entropy_vec
        9/100, # bond_orientation_entropy_vec
        1/100, # coordination_nr_mean_vec
        1/100, # coordination_nr_std_vec 
        ] + [2/100 * 1/13 for _ in range(13) # q_l_mat_values
        ] + [2/100 * 1/13 for _ in range(13) # q_l_mat_uncertainties
        ] + [
        6/100, # vertex_homogeneity_metric_vec
        6/100, # uncoordinated_neighbor_distance_vec
        0.5/100, # ring_size_mean_vec
        0.5/100, # ring_size_std_vec
        1/100, # ring_radius_mean_vec
        1/100, # ring_radius_std_vec
        7/100, # critical_pore_radius_vec
        4.5/100, # anisotropy_metric_from_structure_factor_vec
        4.5/100, # anisotropy_metric_from_structure_factor_bonds_vec
        2/100, # hyperuniformity_alpha_vec
        ])

    return weights


def get_weights_very_focussed():
    """ Set weights for the loss function that focus strongly on bond lengths  
    and angles. The weights are normalized to sum to 1.
    The weights are normalized to sum to 1.
    Metrics describing the primitives: 
    - bond_length_std_vec 30/100
    - bond_angle_std_vec 30/100
    - dihedral_angle_entropy_vec 2/100
    - q_l_mat_values 1/100
    - q_l_mat_uncertainties 1/100
    Metrics describing homogeneity: 
    - vertex_homogeneity_metric_vec 6/100
    - uncoordinated_neighbor_distance_vec 6/100
    - critical_pore_radius_vec 7/100
    Metrics describing anisotropy: 
    - bond_orientation_entropy_vec 6/100
    - anisotropy_metric_from_structure_factor_vec 3/100
    - anisotropy_metric_from_structure_factor_bonds_vec 3/100
    Metrics describing hyperuniformity: 
    - hyperuniformity_alpha_vec 2/100
    Metrics describing ring statistics: 
    - ring_size_mean_vec 0.1/100
    - ring_size_std_vec 0.1/100
    - ring_radius_mean_vec 1/100
    - ring_radius_std_vec 1/100
    Coordination number: 
    - coordination_nr_mean_vec 0.4/100
    - coordination_nr_std_vec 0.4/100
    """

    order_metrics_vec = [
        "bond_length_std_vec", # 30/100
        "bond_angle_std_vec", # 30/100
        "dihedral_angle_entropy_vec", # 2/100
        "bond_orientation_entropy_vec", # 6/100
        "coordination_nr_mean_vec", # 0.4/100
        "coordination_nr_std_vec", # 0.4/100
        "q_l_mat_values", # 1/100 * 1/13
        "q_l_mat_uncertainties", # 1/100 * 1/13
        "vertex_homogeneity_metric_vec", # 6/100
        "uncoordinated_neighbor_distance_vec", # 6/100
        "ring_size_mean_vec", # 0.1/100
        "ring_size_std_vec", # 0.1/100
        "ring_radius_mean_vec", # 1/100
        "ring_radius_std_vec", # 1/100
        "critical_pore_radius_vec", # 7/100
        "anisotropy_metric_from_structure_factor_vec", # 3/100
        "anisotropy_metric_from_structure_factor_bonds_vec", # 3/100
        "hyperuniformity_alpha_vec", # 2/100
        ]

    weights = torch.tensor([
        30/100, # bond_length_std_vec
        30/100, # bond_angle_std_vec
        2/100, # dihedral_angle_entropy_vec
        6/100, # bond_orientation_entropy_vec
        0.4/100, # coordination_nr_mean_vec
        0.4/100, # coordination_nr_std_vec 
        ] + [1/100 * 1/13 for _ in range(13) # q_l_mat_values
        ] + [1/100 * 1/13 for _ in range(13) # q_l_mat_uncertainties
        ] + [
        6/100, # vertex_homogeneity_metric_vec
        6/100, # uncoordinated_neighbor_distance_vec
        0.1/100, # ring_size_mean_vec
        0.1/100, # ring_size_std_vec
        1/100, # ring_radius_mean_vec
        1/100, # ring_radius_std_vec
        7/100, # critical_pore_radius_vec
        3/100, # anisotropy_metric_from_structure_factor_vec
        3/100, # anisotropy_metric_from_structure_factor_bonds_vec
        2/100, # hyperuniformity_alpha_vec
        ])

    return weights


def save_neural_network(neural_network,
             neural_network_data_path="../data/neural_networks/ctn/",
             filename="neural_network.pt"):
    """ Save the neural network to a file. """

    torch.save(neural_network.state_dict(), 
               neural_network_data_path + filename)
    return


def load_neural_network(nr_inputs=3,
             nr_outputs=42,
             nr_neurons=60,
             nr_layers=4,
             nr_neurons_vec=[6,8],
             variable_nr_neurons=False,
             neural_network_data_path="../data/neural_networks/ctn/",
             filename="neural_network.pt"):
    """ Load the neural network from a file. """

    if variable_nr_neurons:
        neural_network = NeuralNetworkVariableNrNeurons(
            nr_inputs=nr_inputs, 
            nr_outputs=nr_outputs,
            nr_neurons_vec=nr_neurons_vec)
    else:
        neural_network = NeuralNetwork(
            nr_inputs=nr_inputs, 
            nr_outputs=nr_outputs, 
            nr_neurons=nr_neurons, 
            nr_layers=nr_layers)
    neural_network.load_state_dict(
        torch.load(neural_network_data_path + filename, weights_only=True))
    neural_network.eval()

    return neural_network


def get_test_loss_and_corr(
        neural_network, 
        weights=None,
        device="cpu", 
        loss_func=None,
        neural_network_data_path="../data/neural_networks/ctn/",
        data_filename="order_metric_dataset.h5",
        nr_outputs=42,
        apply_inverse_pca=False,
        nr_pca_components=10):
    
    if weights is None:
        weights = torch.ones(nr_outputs) / nr_outputs

    if loss_func is None:
        loss_func = WeightedMSE(weights)
    
    # Load data
    trainset, valset, testset = load_dataset(
        neural_network_data_path=neural_network_data_path,
                                          filename=data_filename)
    testloader = DataLoader(testset, shuffle=False)

    loss_sum = 0.0
    data_count = 0

    # Collect predictions and true values for correlation
    all_preds = []
    all_true = []

    with torch.no_grad():
        for data in testloader:
            x_data, y_data = data
            x_data, y_data = x_data.to(device), y_data.to(device)
            outputs = neural_network(x_data)

            if apply_inverse_pca:
                # Apply inverse PCA transformation
                scaler, pca = load_scaler_and_pca(
                    filename=data_filename[:-3]+f"_pca_{nr_pca_components}.h5",
                    neural_network_data_path=neural_network_data_path)
                outputs_np = outputs.detach().cpu().numpy()
                outputs_np = pca.inverse_transform(outputs_np)
                outputs_np = scaler.inverse_transform(outputs_np)
                outputs = torch.tensor(outputs_np, 
                                       dtype=torch.float32).to(device)

            # Accumulate for correlation
            all_preds.append(outputs.cpu().numpy())
            all_true.append(y_data.cpu().numpy())

            # Compute loss
            loss = loss_func(outputs, y_data)
            loss_sum += loss.item()
            data_count += len(x_data)

    test_loss = loss_sum / data_count

    # Concatenate all predictions and true values
    all_preds = np.vstack(all_preds)
    all_true = np.vstack(all_true)

    # Compute Pearson correlation per output dimension
    pearson_correlations = []
    spearman_correlations = []
    r2_scores = []
    for i in range(nr_outputs):
        pearson_corr, _ = pearsonr(all_preds[:, i], all_true[:, i])
        pearson_correlations.append(pearson_corr)
        spearman_corr, _ = spearmanr(all_preds[:, i], all_true[:, i])
        spearman_correlations.append(spearman_corr)
        r2 = r2_score(all_true[:, i], all_preds[:, i])
        r2_scores.append(r2)

    return test_loss, pearson_correlations, spearman_correlations, r2_scores, \
        all_preds, all_true


def train_neural_network(config, 
              weights=None,
              seed=42, 
              checkpoint_path="../data/neural_networks/ctn/checkpoints/",
              neural_network_data_path="../data/neural_networks/ctn/",
              data_filename="order_metric_dataset.h5",
              nr_inputs=8,
              nr_outputs=3,
              variable_nr_neurons=False,
              nr_epochs=10,
              use_checkpoints = False):
    """ Train a neural network with the given configuration. The code is
    adapted from the PyTorch tutorial on hyperparameter tuning:
    https://docs.pytorch.org/tutorials/beginner/hyperparameter_tuning_tutorial.html
    """

    # If the weights are not provided, weigh all output dimensions equally
    if weights is None:
        weights = torch.ones(nr_outputs) / nr_outputs

    np.random.seed(seed)
    torch.manual_seed(seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Define the neural network, loss function, and optimizer
    if variable_nr_neurons:
        neural_network = NeuralNetworkVariableNrNeurons(
            nr_inputs, nr_outputs, config["nr_neurons_vec"]).to(device)
    else:
        neural_network = NeuralNetwork(
            nr_inputs, nr_outputs, config["nr_neurons"], config["nr_layers"]
            ).to(device)
    
    loss_func = WeightedMSE(weights)
    optimizer = optim.SGD(neural_network.parameters(), 
                          lr=config["learning_rate"], 
                          momentum=0.9, 
                          weight_decay=config["weight_decay"])

    # Load from checkpoint if available
    if use_checkpoints:
        checkpoint = get_checkpoint()
        if checkpoint:
            with checkpoint.as_directory() as checkpoint_path:
                checkpoint_data_path = Path(checkpoint_path) / "data.pkl"
                with open(checkpoint_data_path, "rb") as fp:
                    checkpoint_state = pickle.load(fp)
                start_epoch = checkpoint_state["epoch"]
                neural_network.load_state_dict(
                    checkpoint_state["neural_network_state_dict"])
                optimizer.load_state_dict(
                    checkpoint_state["optimizer_state_dict"])
        else:
            start_epoch = 0
    else:
        start_epoch = 0

    # Load data
    trainset, valset, testset = load_dataset(
        neural_network_data_path=neural_network_data_path,
        filename=data_filename)
    trainloader = DataLoader(trainset, batch_size=config["batch_size"],
                             shuffle=True)
    valloader = DataLoader(valset, batch_size=config["batch_size"])

    # Train the neural network
    for epoch in range(start_epoch, nr_epochs):  
        print("Epoch: ", epoch)
        running_loss = 0.0
        epoch_steps = 0
        for i, data in enumerate(trainloader, 0):
            inputs, labels = data
            inputs, labels = inputs.to(device), labels.to(device)

            optimizer.zero_grad()

            # Forward + backward + optimize
            outputs = neural_network(inputs)
            loss = loss_func(outputs, labels)
            loss.backward()
            optimizer.step()

            # Print statistics
            running_loss += loss.item()
            epoch_steps += 1
            if i % 200 == 199:
                print(
                    "[%d, %5d] loss: %.3f"
                    % (epoch + 1, i + 1, running_loss / epoch_steps)
                )
                running_loss = 0.0

        # Validation loss
        val_loss = 0.0
        val_steps = 0
        for i, data in enumerate(valloader, 0):
            with torch.no_grad():
                inputs, labels = data
                inputs, labels = inputs.to(device), labels.to(device)

                outputs = neural_network(inputs)

                loss = loss_func(outputs, labels)
                val_loss += loss.cpu().numpy()
                val_steps += 1

        print("Validation loss: ", val_loss / val_steps)

        # Save checkpoint
        if use_checkpoints:
            checkpoint_data = {
                "epoch": epoch,
                "neural_network_state_dict": neural_network.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
            }
            with tempfile.TemporaryDirectory() as checkpoint_path:
                neural_network_data_path = Path(checkpoint_path) / "data.pkl"
                with open(neural_network_data_path, "wb") as fp:
                    pickle.dump(checkpoint_data, fp)

                checkpoint = Checkpoint.from_directory(checkpoint_path)
                train.report(
                    {"loss": val_loss / val_steps},
                    checkpoint=checkpoint,
                )

    print("Finished Training")

    return neural_network
    

def hyperparam_tuning(weights=None,
        nr_layers_tune=[i for i in range(2, 6)],
        nr_neurons_tune=[i for i in range(40, 100)],
        learning_rate_tune=[0.05, 0.1, 0.2 ],
        batch_size_tune=[1,2,4],
        weight_decay_tune=[0],
        seed=42, 
        checkpoint_path = "../data/neural_networks/ctn/checkpoints/",
        neural_network_data_path = "../data/neural_networks/ctn/",
        data_filename="order_metric_dataset.h5",
        nr_inputs=3,
        nr_outputs=42,
        variable_nr_neurons=False,
        max_nr_epochs=10,
        nr_samples_hyperparams=10,):
    """ Perform hyperparameter tuning using Ray Tune. The code is adapted from 
    the PyTorch tutorial on hyperparameter tuning:
    https://docs.pytorch.org/tutorials/beginner/hyperparameter_tuning_tutorial.html
    """
    
    if weights is None:
        weights = torch.ones(nr_outputs) / nr_outputs
    
    if variable_nr_neurons:
        config = {
            "nr_neurons_vec": tune.choice(nr_neurons_vec_tune),
            "learning_rate": tune.choice(learning_rate_tune),
            "batch_size": tune.choice(batch_size_tune),
            "weight_decay": tune.choice(weight_decay_tune),
        }
    else:
        config = {
            "nr_layers": tune.choice(nr_layers_tune),
            "nr_neurons": tune.choice(nr_neurons_tune),
            "learning_rate": tune.choice(learning_rate_tune),
            "batch_size": tune.choice(batch_size_tune),
            "weight_decay": tune.choice(weight_decay_tune),
        }
    scheduler = ASHAScheduler(
        metric="loss",
        mode="min",
        max_t=max_nr_epochs,
        grace_period=1,
        reduction_factor=2,
    )
    result = tune.run(
        partial(train_neural_network, 
                weights=weights,
                seed=seed, 
                checkpoint_path=checkpoint_path, 
                neural_network_data_path=neural_network_data_path,
                data_filename=data_filename, 
                nr_inputs=nr_inputs, 
                nr_outputs=nr_outputs,
                variable_nr_neurons=variable_nr_neurons,
                use_checkpoints=True),
        config=config,
        num_samples=nr_samples_hyperparams,
        scheduler=scheduler,
    )

    best_trial = result.get_best_trial("loss", "min", "last")
    print(f"Best trial config: {best_trial.config}")
    print(f"Best trial final validation loss: {best_trial.last_result['loss']}")

    if variable_nr_neurons:
        best_trained_model = NeuralNetworkVariableNrNeurons(nr_inputs, nr_outputs, 
                                         config["nr_neurons_vec"])
    else:
        best_trained_model = NeuralNetwork(nr_inputs, nr_outputs, 
                                 best_trial.config["nr_neurons"], 
                                 best_trial.config["nr_layers"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    best_trained_model.to(device)

    best_checkpoint = result.get_best_checkpoint(trial=best_trial, metric="loss", mode="min")
    with best_checkpoint.as_directory() as checkpoint_path:
        network_data_path = Path(checkpoint_path) / "data.pkl"
        with open(network_data_path, "rb") as fp:
            best_checkpoint_data = pickle.load(fp)

        best_trained_model.load_state_dict(best_checkpoint_data["neural_network_state_dict"])
        test_loss, pearson_correlations, spearman_correlations, r2_scores, \
         all_preds, all_true = get_test_loss_and_corr(
            best_trained_model, 
            weights=weights,
            device=device, 
            neural_network_data_path=neural_network_data_path,
            data_filename=data_filename,
            nr_outputs=nr_outputs)
        print("Best trial test set loss: {}".format(test_loss))
