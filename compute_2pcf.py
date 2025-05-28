# =========================================================================== #
"""
Created on May 27, 2025 
@authors: Zachery Brown

Description:
    
"""
# =========================================================================== #

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
import json
import argparse
import src.utils as utils
import src.pair_counter as pc

# =========================================================================== #


def main():
    parser = argparse.ArgumentParser(
        description="Compute the 2pcf with parameters specified "+\
                    "by a configuration file.")
    parser.add_argument(
        "-c", "--config",
        type=str,
        required=True,
        help="Path to the configuration file"
    )

    args = parser.parse_args()
    with open(args.config) as cfgDict:
        cfg_set = json.load(cfgDict)
        
    print(f"Computing 2pcf for file {cfg_set['data_raw_file']}")
    print(f"Randoms file {cfg_set['rand_raw_file']}")
    
    # Load the data coords
    data_raw = fits.open(cfg_set['data_raw_file'])[1].data

    # Load the random coords
    rand_raw = fits.open(cfg_set['rand_raw_file'])[1].data

    # Set the z-limits
    zlims = (cfg_set['zlims'][0],cfg_set['zlims'][1])
    
    print('Pruning catalogs...')
    
    # Prune the data catalog
    data = utils.catalog_slicer(
        data_raw,zlims,resample_factor=cfg_set['data_resample_factor'])
    
    print(f'Original data catalog had {len(data_raw)} objects')
    print(f'Sliced data catalog has {len(data)} objects')
    print('Data catalog is '+\
          f'{np.round(100*(len(data)/len(data_raw)),decimals=2)}'+\
          ' % of original')

    # Prune the randoms catalog
    rand = utils.catalog_slicer(
        rand_raw,zlims,resample_factor=cfg_set['rand_resample_factor'])

    print(f'Original rand cat had {len(rand_raw)} objects')
    print(f'Sliced rand cat has {len(rand)} objects')
    print('Rand catalog is '+\
          f'{np.round(100*(len(rand)/len(rand_raw)),decimals=2)}'+\
          ' % of original')

    print('Prepping binning...')
    
    # Specify cosmology
    cosmo = FlatLambdaCDM(H0=cfg_set['H0'], Om0=cfg_set['Om0'])

    # Number of threads to use
    nthreads = cfg_set['nthreads']

    # Create the bins array
    rbins = np.linspace(cfg_set['rbins'][0],
                        cfg_set['rbins'][1],
                        cfg_set['rbins'][2])

    # Specify the max. of the cosine of the angle to the LOS
    # for DD(s, mu)
    mu_max = cfg_set['mu_max']

    # Specify the number of linear bins in `mu`
    Nmu_bins = cfg_set['Nmu_bins']

    # Set additional arguments to pass to corrfunc's DDsmu_mocks
    corrfunc_kwargs = cfg_set['corrfunc_kwargs']
    
    # Count pairs
    print('Counting DD pairs...')
    DD = pc.count_DD_pairs(data,cosmo,nthreads,rbins,mu_max,Nmu_bins,
                           corrfunc_kwargs)

    print('Counting DR pairs...')
    DR = pc.count_DR_pairs(data,rand,cosmo,nthreads,rbins,mu_max,Nmu_bins,
                           corrfunc_kwargs)

    print('Counting RR pairs...')
    RR = pc.count_RR_pairs(rand,cosmo,nthreads,rbins,mu_max,Nmu_bins,
                           corrfunc_kwargs)

    print('Wrapping results...')

    # Set desired ells
    ell = cfg_set['ell']

    s, xi_ell = utils.corr_2pcf_legendre(
        data,rand,ell,mu_max,Nmu_bins,DD,DR,RR)

    utils.write_2pcf_to_file(s, xi_ell, cfg_set['corr_filename'])

    utils.plot_2pcf(s, xi_ell, cfg_set['plot_filename'])

if __name__ == "__main__":
    main()

