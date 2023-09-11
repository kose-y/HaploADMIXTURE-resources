# HaploADMIXTURE-resources

This repository contains information for reproducing experiments in the HaploADMIXTURE manuscript:

_To be updated_

To download the 1000 Genomes Project (TGP), Human Genome Diversity Project (HGDP), and Human Origins (HO) public data sets, follow the steps in https://github.com/sriramlab/SCOPE/tree/master/misc/real_data and https://github.com/StoreyLab/terastructure/tree/master/scripts/data. 

The directory `environments` defines Julia package environment for the experiments. To install all the packages, 
- launch Julia with `--project=<path-to-this directory>/environments/v1.7` or activate the environment with `using Pkg; Pkg.activate(normpath("<path-to-this directory>/environments/v1.7"))`,
- then run `using Pkg; Pkg.instantiate()`.
- Note that the experiments were conducted in Julia 1.7, while the Manifest information in the published packages has been upgraded to Julia v1.9.
