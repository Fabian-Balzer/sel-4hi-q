# Sel-4Hi-Q
The code supplied here is capable of processing the data provided by the different catalogues and to perform the photometric redshift fitting via LePhare.
In addition to that, various plotting routines are provided to perform analysis of the joint input catalogue and the output coming from LePhare.


## Requirements
The following programs are necessary for the code to work:

- The current version of [LePhare++](https://gitlab.lam.fr/Galaxies/LEPHARE), along with proper setup of the `LEPHAREDIR` and `LEPHAREWORK` environment variables.
- A python 3.8+ installation, including pandas, matplotlib and numpy (available i. e. via Anaconda).
- A directory containing the input catalogues in a directory with the environment name `CATPATH`.
- A directory you wish to work in with the environment name `LEPHARE`.
After the repository has been cloned, the program `startup.py` should be run using

`>>> python startup.py`

This routine sets up the proper directory structure, provides information about missing components, and it also resets the basic configuration files `general_config.ini` (which contains general information about the system and should almost never be altered) and `custom_config.ini`.
In addition to that, it downloads the program [jystilts](http://www.star.bris.ac.uk/~mbt/stilts/), which is the \codep{jython} installation of STILTS that is used to process much of the data.
Note: After additional components have been installed and new environment variables have been set, the `startup.py`  program needs to be run again!

## Running the code
The whole code can then be run from the main directory using

`>>> python run_all.py`

This will execute the code based on the parameters set in `custom_config.ini` (please refer to the thesis for an explanation of the parameters).
If multiple different configuration files are desired to be used, they need to be stored in the `config` directory and the active one needs to be specified in `general_config.ini`, changing the parameter for `current_config` to the name of the new configuration file.