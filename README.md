# shotgun-pipeline
PennCHOP shotgun metagenomics pipeline implemented in snakem

### snakemake setup
Do this once:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda create --name shotgun --file requirements.txt

source activate shotgun-pipeline

# conda only keep track of the packages it installed. TODO: create conda package
pip install -i https://pypi.anaconda.org/zhaoc1/simple dnabc
pip install -i https://pypi.anaconda.org/zhaoc1/simple illqc
pip install -i https://pypi.anaconda.org/zhaoc1/simple decontam

snakemake
```

### Development

To updated the requirements file (after installing some new package):
```
conda list --name shotgun --explicit > requirements.txt
```

To update your conda environment with a new requirements file:
```
conda install --name shotgun --file requirements.txt
```