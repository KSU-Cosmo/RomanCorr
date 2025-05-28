"""
Created on Thu May 22 11:54:26 2025

@author: Zachery Brown

"""

# =========================================================================== #

import numpy as np
from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
import astropy.units as u
import astropy

# =========================================================================== #

def count_DD_pairs(data: astropy.io.fits.fitsrec.FITS_rec, 
                   cosmology: astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM,
                   nthreads: int,
                   rbins: np.ndarray,
                   mu_max: float,
                   Nmu_bins: int,
                   corrfunc_kwargs: dict):
    """
    !=!=!=!=!=!=!=!=!=!=!=!=!=!=! TO DO !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    --> Include non-flat, non Lambda CDM functionality
    
    Using corrfunc, count the number of 'data'/'data' pairs in the selected
    catalog.
    
    Parameters
    -----------
    
    data: astropy.io.fits.fitsrec.FITS_rec, required
        Catalog with column format ['ra','dec','z','wts'].

    cosmology: astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM, required
        A cosmology from astropy's cosmology class.

    nthreads: int, required
        The number of threads to pass to corrfunc's pair counter engine.

    rbins: numpy.ndarray, required
        The bin edges in r (or s) used to compute the 2pcf.

    mu_max: float, required
        The largest value of mu used to compute the 2pcf. mu is defined 
        to be the cosine of the angle with respect to the line of sight.

    Nmu_bins: int, required
        The number of bins in mu, begining from 0 and going to mu_max.

    corrfunc_kwargs: dict, required
        Additional arguments to pass to corrfunc's DDsmu_mocks. If no 
        other arguments are desired, the dictionary can be left empty.
        For more details regarding corrfunc options, please see their
        API reference at: https://corrfunc.readthedocs.io
        
    Returns
    --------

    DD_pairs: numpy.ndarray
        The raw output of corrfunc's pair counting engine. Array contains
        the number of pairs at a given r and mu.

    """
    
    # Pass arguments to DDsmu_mocks and count pairs
    DD_pairs = DDsmu_mocks(
        autocorr=1,
        cosmology=1,
        nthreads=nthreads,
        mu_max=mu_max,
        nmu_bins=Nmu_bins,
        binfile=rbins,
        RA1=np.asarray(data['ra'],dtype=np.float64),
        DEC1=np.asarray(data['dec'],dtype=np.float64),
        CZ1=np.asarray(cosmology.comoving_distance(
            data['z']).to(u.Mpc).value*cosmology.h,dtype=np.float64),
        weights1=np.asarray(data['wts'],dtype=np.float64),
        is_comoving_dist=True,
        output_savg=True,
        **corrfunc_kwargs)

    # Return the corrfunc output array
    return DD_pairs


def count_RR_pairs(rand: astropy.io.fits.fitsrec.FITS_rec, 
                   cosmology: astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM,
                   nthreads: int,
                   rbins: np.ndarray,
                   mu_max: float,
                   Nmu_bins: int,
                   corrfunc_kwargs: dict):
    """
    !=!=!=!=!=!=!=!=!=!=!=!=!=!=! TO DO !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    --> Include non-flat, non Lambda CDM functionality
    
    Using corrfunc, count the number of 'random'/'random' pairs in the 
    selected catalog.

    Note: This function is identical to count_DD_pairs. I have chosen to
          redefine it as a new function for clarity, so the user is not 
          using a function with 'DD' in the name to count random-random
          pairs.
          
    Parameters
    -----------
    
    rand: astropy.io.fits.fitsrec.FITS_rec, required
        Catalog with column format ['ra','dec','z','wts'].

    cosmology: astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM, required
        A cosmology from astropy's cosmology class.

    nthreads: int, required
        The number of threads to pass to corrfunc's pair counter engine.

    rbins: numpy.ndarray, required
        The bin edges in r (or s) used to compute the 2pcf.

    mu_max: float, required
        The largest value of mu used to compute the 2pcf. mu is defined 
        to be the cosine of the angle with respect to the line of sight.

    Nmu_bins: int, required
        The number of bins in mu, begining from 0 and going to mu_max.

    corrfunc_kwargs: dict, required
        Additional arguments to pass to corrfunc's DDsmu_mocks. If no 
        other arguments are desired, the dictionary can be left empty.
        For more details regarding corrfunc options, please see their
        API reference at: https://corrfunc.readthedocs.io          
        
    Returns
    --------

    RR_pairs: numpy.ndarray
        The raw output of corrfunc's pair counting engine. Array contains
        the number of pairs at a given r and mu.

    """
    
    # Pass arguments to DDsmu_mocks and count pairs
    RR_pairs = DDsmu_mocks(
        autocorr=1,
        cosmology=1,
        nthreads=nthreads,
        mu_max=mu_max,
        nmu_bins=Nmu_bins,
        binfile=rbins,
        RA1=np.asarray(rand['ra'],dtype=np.float64),
        DEC1=np.asarray(rand['dec'],dtype=np.float64),
        CZ1=np.asarray(cosmology.comoving_distance(
            rand['z']).to(u.Mpc).value*cosmology.h,dtype=np.float64),
        weights1=np.asarray(rand['wts'],dtype=np.float64),
        is_comoving_dist=True,
        output_savg=True,
        **corrfunc_kwargs)
    
    # Return the corrfunc output array
    return RR_pairs


def count_DR_pairs(data: astropy.io.fits.fitsrec.FITS_rec,
                   rand: astropy.io.fits.fitsrec.FITS_rec,
                   cosmology: astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM,
                   nthreads: int,
                   rbins: np.ndarray,
                   mu_max: float,
                   Nmu_bins: int,
                   corrfunc_kwargs: dict):
    """
    !=!=!=!=!=!=!=!=!=!=!=!=!=!=! TO DO !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    --> Include non-flat, non Lambda CDM functionality
    
    Using corrfunc, count the number of 'data'/'random' pairs in the 
    selected catalogs.
          
    Parameters
    -----------
    
    data: astropy.io.fits.fitsrec.FITS_rec, required
        Catalog with column format ['ra','dec','z','wts'].
    
    rand: astropy.io.fits.fitsrec.FITS_rec, required
        Catalog with column format ['ra','dec','z','wts'].

    cosmology: astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM, required
        A cosmology from astropy's cosmology class.

    nthreads: int, required
        The number of threads to pass to corrfunc's pair counter engine.

    rbins: numpy.ndarray, required
        The bin edges in r (or s) used to compute the 2pcf.

    mu_max: float, required
        The largest value of mu used to compute the 2pcf. mu is defined 
        to be the cosine of the angle with respect to the line of sight.

    Nmu_bins: int, required
        The number of bins in mu, begining from 0 and going to mu_max.

    corrfunc_kwargs: dict, required
        Additional arguments to pass to corrfunc's DDsmu_mocks. If no 
        other arguments are desired, the dictionary can be left empty.
        For more details regarding corrfunc options, please see their
        API reference at: https://corrfunc.readthedocs.io          
        
    Returns
    --------

    DR_pairs: numpy.ndarray
        The raw output of corrfunc's pair counting engine. Array contains
        the number of pairs at a given r and mu.

    """
    
    # Pass arguments to DDsmu_mocks and count pairs
    DR_pairs = DDsmu_mocks(
        autocorr=0,
        cosmology=1,
        nthreads=nthreads,
        mu_max=mu_max,
        nmu_bins=Nmu_bins,
        binfile=rbins,
        RA1=np.asarray(data['ra'],dtype=np.float64),
        DEC1=np.asarray(data['dec'],dtype=np.float64),
        CZ1=np.asarray(cosmology.comoving_distance(
            data['z']).to(u.Mpc).value*cosmology.h,dtype=np.float64),
        weights1=np.asarray(data['wts'],dtype=np.float64),
        RA2=np.asarray(rand['ra'],dtype=np.float64),
        DEC2=np.asarray(rand['dec'],dtype=np.float64),
        CZ2=np.asarray(cosmology.comoving_distance(
            rand['z']).to(u.Mpc).value*cosmology.h,dtype=np.float64),
        weights2=np.asarray(rand['wts'],dtype=np.float64),
        is_comoving_dist=True,
        output_savg=True,
        **corrfunc_kwargs)
    
    # Return the corrfunc output array
    return DR_pairs

# =========================================================================== #

