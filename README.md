# RomanCorr
Some wrapper functions to compute clustering statistics for Roman galaxy survey data.

## Installation Instructions

### Navigate to desired directory
`cd /path/to/your/directory`

### Name the venv
`python3 -m venv RomanCorr-env`

### Activate the venv
`source RomanCorr-env/bin/activate`

### Upgrade the pip install wheel
`pip install --upgrade pip setuptools wheel`

### Install the requirements
`pip install -r requirements.txt`

### Install corrfunc dependencies seperately
`pip install cython setuptools`

`pip install corrfunc`

## Computing the 2pcf

`python compute_2pcf.py --config ./configs/sample_config.txt`

## Using jupyter in this virtual environment

### Name and register the kernel for notebooks
`python -m ipykernel install --user --name=RomanCorr --display-name "Python (RomanCorr)"`

### Lauch jupyter and open desired notebook
`jupyter notebook`
