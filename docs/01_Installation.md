# Installation
## Install and configure conda package manager
### Install Miniconda (optional)
We recommend using conda to install STELR's software dependencies. If your system doesn't have conda installed, please use the following steps to install Miniconda (Python 3.X).
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc

conda init # this step requires you to close and open a new terminal before it take effect
conda update conda # update conda

### Install mamba (optional)
Mamba is a reimplementation of the conda package manager in C++. There is significant speed improvement on STELR installation using mamba versus conda. Please use following command to install mamba into the base conda environment. You can skip this step if you already have mamba installed.
```
conda install mamba -n base -c conda-forge
```
For more on mamba: see [Mamba's documentation](https://mamba.readthedocs.io/en/latest/).

STELR and all its software dependencies can be installed directly using the STELR git repository. Note: installation using this approach ensures fixed dependency versions.
```
git clone git@github.com:bergmanlab/STELR.git
cd STELR
mamba env create -f envs/stelr.yaml
conda activate STELR
```
## Activate TELR Conda Environment
The STELR conda environment must always be activated prior to running STELR. This step adds STELR's dependencies installed in the STELR conda environment to the environment PATH.
```
conda activate STELR
```
NOTE: Sometimes activating conda environments does not work via conda activate env when run through a script submitted to a queueing system, this can be fixed by activating the environment in the script as shown below.
```
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate STELR
```
For more on Conda: see the [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/index.html).

## Run STELR on test dataset
A test dataset is provided in the `test/` directory, you can test whether your STELR installation is successful by cloning STELR repository and running STELR on the test dataset within the local STELR repository. The test run should generally take less than one minute to finish.
```
git clone git@github.com:bergmanlab/STELR.git
cd STELR/test
conda activate STELR
python3 ../src/telr/stelr.py -o test_output -i reads.fasta -r ref_38kb.fasta -l library.fasta
```