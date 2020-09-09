# Installation
TELR is written in python3 and is designed to run on linux operating systems. Installation of software dependencies for McClintock is automated by Conda, thus a working installation of Conda is required to install McClintock. Conda can be installed via the Miniconda installer.

### Installing Miniconda (Python 3.X)
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
```

### Update Conda
```
conda update conda
```

## Install dependencies
After installing and updating Conda, you can now use conda to install dependencies and create running environment for TELR.
```
conda create -n TELR_env python=3.6 -y
conda activate TELR_env
conda install -c conda-forge -y biopython
conda install -y pandas=1.0.0
conda install -y repeatmasker=4.0.7
conda install -y samtools=1.9
conda install -y bcftools=1.9
conda install -y bedtools
conda install -y ngmlr=0.2.7
conda install -y sniffles=1.0.12
conda install -y wtdbg
conda install -y seqtk
conda install -y minimap2
conda install -y svim=1.3
conda install -y mafft
conda install -y raxml
conda install -c bioconda iqtree
```

## Quick Start
```
git clone git@github.com:bergmanlab/TELR.git
cd TELR
./telr.py --help
```

## Running on a test dataset
```
conda activate TELR_env
cd test
../telr.py -o test_output -i reads.fasta -r ref_38kb.fasta -l library.fasta
```