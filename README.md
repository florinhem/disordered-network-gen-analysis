This repository contains code for generating and analyzing spatial networks with arbitrary coordination statistics using an extended Wooten–Weaire–Winer (WWW) algorithm. The theoretical background, equations, and applications are described in:

> [**Algorithmic Design of Disordered Networks With Arbitrary Coordination: Application to Biophotonics**  
> *Advanced Functional Materials*](https://advanced.onlinelibrary.wiley.com/doi/10.1002/adfm.202600037)

## Repository Structure

- **`julia_code/`**  
  Contains the core Julia implementation for generating and analyzing spatial networks.  
  An overview of the functionality is provided in the notebook  
  **`networks_walkthrough.ipynb`**.

  The datasets presented in the article were generated using the script  
  **`create_dataset.jl`**.  
  For practical reasons, the number of network samples in this script is reduced compared to the full datasets used in the publication.

- **`python_order_metric_prediction/`**  
  Contains the Python code for predicting order metrics of networks generated from the initial __ctn__ network.  
  The workflow is demonstrated in the notebook  
  **`order_metric_prediction.ipynb`**.
