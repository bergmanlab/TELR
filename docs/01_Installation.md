# Installation
## Use Conda to install software dependencies
TELR is written in python3 and is designed to run on linux operating system. Installation of software dependencies for TELR is automated by Conda, thus a working installation of Conda is required to install TELR. Conda can be installed via the Miniconda installer.

### Installing Miniconda (Python 3.X)
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
conda init
```
- `conda init` requires you to close and open a new terminal before it take effect

### Update Conda
```
conda update conda
```

## Install software dependencies
After installing and updating Conda, you can now use conda to install dependencies and create running environment for TELR.
### Clone TELR Repository
```
git clone git@github.com:bergmanlab/TELR.git
cd TELR
```
### Create TELR Conda Environment
```
conda env create -f envs/telr.yml
```
- This installs all the software dependencies needed to run TELR into the TELR Conda environment.

### Activate TELR Conda Environment
```
conda activate TELR_env
```
- This adds the dependencies installed in the TELR conda environment to the environment PATH so that they can be used by the TELR scripts.
- This environment must always be activated prior to running any of the TELR scripts
- NOTE: Sometimes activating conda environments does not work via conda activate myenv when run through a script submitted to a queueing system, this can be fixed by activating the environment in the script as shown below
```
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate TELR_env
```
- For more on Conda: see the [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/index.html).

## Running TELR on test dataset
- A test dataset is provided in the `test/` directory, you can test whether your TELR installation is successful by running TELR on this dataset, which should take less than one minute to finish on a single thread machine.
```
conda activate TELR_env
cd test
../telr.py -o test_output -i reads.fasta -r ref_38kb.fasta -l library.fasta
```