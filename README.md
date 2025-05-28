# RomanCorr
Some wrapper functions to compute clustering statistics for Roman galaxy survey data.

## Installation Instructions

The instructions here create a python virtual environment called RomanCorr-env. I have found that it is easier to install corrfunc and its dependencies seperately.

Note: corrfunc requires a local c compiler!

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

Below is how to compute the 2pcf. You may edit the configuration file and make as many as you'd like. When running this routine, simply specify the path to your config file from the command line.

You will want to set the desired paramters such as filenames and binning in the configuration file. Please note, in the config. file there are two parameters called `data_resample_factor` and `rand_resample_factor`. These allow you to randomly resample a catalog to a smaller size for fast computation and diagnostics. Please note that in the sample configuration file these are both set to 1%. You may wish to change this!

Note: There are no data catalogs provided in this repo. They are simply too large. These routines will work for any catalog, provided it is stored as an astropy .fits file. It MUST have the following column names `ra` (right ascension), `dec` (declination), `z` (redshift), and `wts` (weights). If you would like to run without weights, simply populate the `wts` column with all ones.

`python compute_2pcf.py --config ./configs/sample_config.txt`

## Using jupyter in this virtual environment

If you would prefer to use jupyter notebooks, the following lines will register a kernel in the RomanCorr-env environment, allowing you to run these routines interactively, and do development work. 

### Name and register the kernel for notebooks
`python -m ipykernel install --user --name=RomanCorr --display-name "Python (RomanCorr)"`

### Lauch jupyter and open desired notebook
`jupyter notebook`
