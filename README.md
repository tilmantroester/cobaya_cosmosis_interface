# An interface to run CosmoSIS likelihoods in cobaya

This package provides a cobaya `Likelihood` that reads in a CosmoSIS pipeline from a `.ini` file. 
It provides two modes that differ in where the Boltzmann solver runs:
- Run the CosmoSIS pipeline as-is and ignore any cobaya theory code that might also be running.
- Remove the Boltzmann modules (e.g. CAMB) from the CosmoSIS pipeline and get the necessary information from the cobaya theory code.

The first option reproduces the CosmoSIS log-likelihoods exactly, while the 2nd option makes more sense when combining with other cobaya likelihoods. 
This is controlled with the `use_cobaya_theory` option in the likelihood class.

The package also provides some utility functions in `cosmosis_cobaya_interface.prepare_config` that take the CosmoSIS configurations and create the corresponding cobaya `.yaml` file. 
In particular, it reads the `values` and `prior` files and tries to create the corresponding cobaya priors. 


