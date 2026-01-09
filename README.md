Code to generate and analyze spatial networks with the extended Wooten-Weaire-Winer algorithm.

The `julia_code` directory contains the main code for generating and analyzing spatial networks. The features of this code are explained in the notebook `networks_walkthrough.ipynb`. The datasets in the article were generated using the script `create_dataset.jl`. To reduce the runtime of this script, the number of network samples is chosen to be smaller than in the article dataset.

In the `python_order_metric_prediction` directory, you find the code to predict order metrics for networks generated from the initial \textbf{ctn} network. Here, the workflow is demonstrated in the notebook `order_metric_prediction.ipynb`.
