# HaploADMIXTURE-resources

This repository contains information for reproducing experiments in the HaploADMIXTURE manuscript:

_To be updated_

To download the 1000 Genomes Project (TGP), Human Genome Diversity Project (HGDP), and Human Origins (HO) public data sets, follow the steps in https://github.com/sriramlab/SCOPE/tree/master/misc/real_data and https://github.com/StoreyLab/terastructure/tree/master/scripts/data. 

The files `Project.toml` and `Manifest.toml` define Julia package environment for the experiments. To install all the packages, 
- launch Julia with `--project=<path-to-this directory>` or activate the environment with `using Pkg; Pkg.activate(normpath("<path-to-this directory>"))`,
- then run `using Pkg; Pkg.instantiate()`.

_Please note that the name of the SKFR package has changed from SKFR to SparseKmeansFeatureRanking._


Each real data directory (TGP, HGDP, HO) contains the following:
- `filter_[].jl`: filtering the original datasets using SKFR.
- `run_CV.jl`: (4-fold) cross-validation for selecting the value of `K`, the number of populations, the training step.
- `run_CV_test.jl`: (4-fold) cross-validation for selecting the value of `K`, the number of populations, the test step.
- `run_[]_Ks.jl`: Run the HaploADMIXTURE model for different values of `K`, to select `K` using Akaike information criterion.
- `run_[]_adm.jl`: Run the original admixture model using OpenADMIXTURE.
- `run_[]_full.jl`: Run HaploADMIXTURE with all the SNPs in the original dataset.
- `run_[].jl`: Run HaploADMIXTURE with the dataset filtered with `filter_[].jl`.
- `run_admix_full.jl`: Run OpenADMIXTURE with the original dataset.

Please note that the released version of the package on https://github.com/OpenMendel/HaploADMIXTURE.jl contains the function `run_admixture()` to simplify the code.

For the `Simulation` directory, there is the `simulation.jl` file to generate datasets for simulation.
