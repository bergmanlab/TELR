# Installation
## Install Miniconda (optional)
We recommened using conda to install TELR and its software dependencies. If your system doesn't have conda installed, please use the following steps to install Miniconda (Python 3.X).
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc

conda init # this step requires you to close and open a new terminal before it take effect
conda update conda # update conda
```
## Set up conda channels (optional)
TELR is hosted under bioconda channel (https://anaconda.org/bioconda/telr). After installing conda you will need to add the bioconda channel as well as the other channels bioconda depends on. You can skip this step if you already have conda installed and bioconda channel configured.
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
## Install TELR
TELR and all its software dependencies can be installed using conda. We recommend installing TELR in a new conda environment.
```
conda create -n TELR --channel bioconda telr
```
## Activate TELR Conda Environment
The TELR conda environment must always be activated prior to running TELR. This step adds TELR and its dependencies installed in the TELR conda environment to the environment PATH so that they can be used by TELR.
```
conda activate TELR
```
NOTE: Sometimes activating conda environments does not work via conda activate env when run through a script submitted to a queueing system, this can be fixed by activating the environment in the script as shown below.
```
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate TELR
```
For more on Conda: see the [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/index.html).

## Run TELR on test dataset
A test dataset is provided in the `test/` directory, you can test whether your TELR installation is successful by cloning TELR repository and running TELR on the test dataset within local TELR repository. The test run should generally take less than one minute to finish.
```
git clone git@github.com:bergmanlab/TELR.git
cd TELR/test
conda activate TELR
telr -o test_output -i reads.fasta -r ref_38kb.fasta -l library.fasta
```