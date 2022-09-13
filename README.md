# synthetic-key-steps
Code repository for "Computer-Aided Key Step Generation in Alkaloid Total Synthesis"
Contains notebooks for generation and plotting of graph edit distances, bond type evaluation, visualization of all connectivity matrices, and plots of other molecular fingerprint-based similarity metrics. 
A full printout of package versions in the conda environment used can be found in software_versions.txt. Select version numbers that may be of particular interest are: ipython 7.27.0, jupyterlab 3.1.12, matplotlib 3.4.3, numpy 1.20.3, pandas 1.3.3, RDKit 2021.03.5. These need not be strictly adhered to.

A guide on RDKit installation can be found here: https://www.rdkit.org/docs/Install.html
Notebooks beginning with 01 build matrix-encoded synthetic routes, given a target adjacency matrix, and a .csv file detailing each step’s atom additions and bond edits. 
Notebooks beginning with 02 evaluate and plot bond edit distances from connectivity matrices.
03_plot_matrices.ipynb evaluates bond types for a given route and plots individual connectivity matrices.
04_other_simi_metrics.ipynb computes and compares various molecular fingerprint based similarity metrics, using selected stemoamide synthesis routes.

To run the notebooks presented here, create a conda environment with the appropriate package versions (the RDKit installation guide linked above has a walkthrough), then download this repository into a suitable location on the local machine.
