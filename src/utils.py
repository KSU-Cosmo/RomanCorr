"""
Created on Thu May 22 11:54:26 2025

@author: Zachery Brown

"""

# =========================================================================== #

import astropy
import numpy as np
import numpy
from numpy.polynomial.legendre import legval
from matplotlib import pyplot as plt
from astropy.io import fits

# =========================================================================== #


def catalog_slicer(catalog: astropy.io.fits.fitsrec.FITS_rec, 
                   zlims: tuple, resample_factor: float = None):
    """
    Select objects from a fits catalog in a desired redshift range, and 
    randomly resample objects to reduce catalog size if desired.

    Parameters
    -----------

    catalog: astropy.io.fits.fitsrec.FITS_rec, required
        Standard catalog format is columns ['ra','dec','z','wts'], 
        however this function requires only the 'z' column.

    zlims: tuple, required
        A tuple of the limits of the redshift selection. Redshift limits
        should be in ascending order!

    resample_factor, float, optional
        The fraction of objects from the catalog to keep in the selection.
        If None (default), all objects within the redshift range will be
        passed to the output. Other wise, objects will be randomly resampled.

    Returns
    --------

    zsliced_catalog: astropy.io.fits.fitsrec.FITS_rec
        Returns a catalog with the same columns as the inputs, keeping
        only the objects that fall within the desired redshift range.

    """
    
    # If there is no random resampling
    if resample_factor is None:
        
        # Perform z-slicing
        z_selection = (catalog['z']>zlims[0])&\
                      (catalog['z']<zlims[1])
        
        # Return new catalog with z-selection
        return catalog[z_selection]

    # If random resampling is desired
    elif resample_factor is not None:

        # Size of the resampled catalog
        catalog_subsampled_N = int(resample_factor*len(catalog))

        # Indices of selected objects
        catalog_subsampled_ids = np.random.choice(
            np.arange(len(catalog), dtype=int),
            size=catalog_subsampled_N,
            replace=False)

        # Reduce catalog size from resampling
        catalog_resampled = catalog[catalog_subsampled_ids]
        
        # Perform z-slicing
        z_selection = (catalog_resampled['z']>zlims[0])&\
                      (catalog_resampled['z']<zlims[1])

        # Return new catalog with z-selection and resampling
        return catalog_resampled[z_selection]
    

def corr_2pcf_legendre(data: astropy.io.fits.fitsrec.FITS_rec,
                       rand: astropy.io.fits.fitsrec.FITS_rec,
                       ell: list, mu_max: float, Nmu_bins: int, 
                       DD_pairs: numpy.ndarray, DR_pairs: numpy.ndarray, 
                       RR_pairs: numpy.ndarray):
    """
    Compute the LS estimator of the 2pcf, and perform the Legendre
    decomposition, returning the desired multipoles.
    
    Parameters
    -----------

    data: astropy.io.fits.fitsrec.FITS_rec, required
        Catalog with column format ['ra','dec','z','wts']. For this
        routine only the weights column is required for normalization.
    
    rand: astropy.io.fits.fitsrec.FITS_rec, required
        Catalog with column format ['ra','dec','z','wts']. For this
        routine only the weights column is required for normalization.

    ell: list, required
        A list of the desired multipoles of the 2pcf to return. Generally
        only the even multipoles are meaningful and non-zero.

    mu_max: float, required
        The largest value of mu used to compute the 2pcf. mu is defined 
        to be the cosine of the angle with respect to the line of sight.

    Nmu_bins: int, required
        The number of bins in mu, begining from 0 and going to mu_max.

    DD_pairs: numpy.ndarray, required
        The raw output of corrfunc's pair counting engine. Array contains
        the number of pairs at a given r and mu.

    DR_pairs: numpy.ndarray, required
        The raw output of corrfunc's pair counting engine. Array contains
        the number of pairs at a given r and mu.
        
    RR_pairs: numpy.ndarray, required
        The raw output of corrfunc's pair counting engine. Array contains
        the number of pairs at a given r and mu.

    Returns
    --------

    s: numpy.ndarray
        The value of s in each bin. Since 'output_savg' is passed to 
        corrfunc's pair counter as True, it will return the precise 
        s value in each bin. Since we are integrating over mu we will
        average the true s values for the output.

    xi_ell: dict
        The 2pcf for each of the desired ells. Individual correlations 
        are stored as arrays.

    """

    # Find the weighted number of objects in data and randoms
    nD=np.sum(data['wts'])
    nR=np.sum(rand['wts'])

    # Define the normalization factor for the pair counts
    norm_DD = 1.0 / (nD **2)
    norm_DR = 1.0 / (nD * nR)
    norm_RR = 1.0 / (nR **2)

    # Set the limits for mu integration and the s-binning
    sbin_lims = np.unique(
        [(row['smin'], row['smax']) for row in DD_pairs], axis=0)
    mu_edges = np.linspace(0, mu_max, Nmu_bins + 1)
    mu_centers = 0.5 * (mu_edges[:-1] + mu_edges[1:])
    dmu = mu_max / Nmu_bins

    # Initialize lists for s-bins and xi_ell
    xi_ell = {l: [] for l in ell}
    s_vals = []

    # Loop over each s-bin
    for si in range(len(sbin_lims)):

        # Mask only the current s-bin and normalize pair counts
        mask = (DD_pairs['smin'] == sbin_lims[si][0])&\
               (DD_pairs['smax'] == sbin_lims[si][1])
        dd = DD_pairs['npairs'][mask].astype(float) * norm_DD
        dr = DR_pairs['npairs'][mask].astype(float) * norm_DR
        rr = RR_pairs['npairs'][mask].astype(float) * norm_RR
        
        # Return the average true s value
        s_avg = np.mean(DD_pairs['savg'][mask])
        s_vals.append(s_avg)

        # Compute the LS estimator of xi
        xi_mu = (dd - 2 * dr + rr) / rr

        # For each l in desired ells
        for l in ell:

            # Integrate with respect to mu
            coeffs = np.zeros(l + 1)
            coeffs[l] = 1
            L_l_mu = legval(mu_centers, coeffs)
            xi_l = np.sum(xi_mu * L_l_mu) * dmu
            xi_ell[l].append(xi_l)

    # Make output into np arrays
    for l in ell:
        xi_ell[l] = np.asarray(xi_ell[l])
    s = np.asarray(s_vals)

    # Return s bins and xi_ell values
    return s, xi_ell
    

def plot_2pcf(s: numpy.ndarray, xi_ell: dict, filename: str = None):
    
    """
    Plot the 2pcf multipoles and save the figure if desired.
    
    Parameters
    -----------

    s: numpy.ndarray, required
        The value of s in each bin of the 2pcf.
    
    xi_ell: dict, required
        The 2pcf for each of the desired ells. Individual correlations 
        are stored as arrays.

    filename: str, optional
        If filename is set, it will attempt to save the figure. Please specify
        the full path to the desired pdf file.
        
    Returns
    --------

    None

    """
    
    # Make axes and plot xi
    fig, ax = plt.subplots(figsize=(5,4))
    for ell in xi_ell.keys():
        ax.plot(s, s**2*xi_ell[ell], label=f'l={ell}',lw=2.4,alpha=.8)
    ax.set(xlabel='s [Mpc/h]',ylabel=r'$\xi_\ell(s)$')
    ax.legend()
    
    # If filename is defined
    if filename is not None:
        plt.savefig(filename)
        
    # Display figure    
    plt.show()
    
    # No return
    return None


def write_2pcf_to_file(s: numpy.ndarray, xi_ell: dict,
                       filename: str):
    
    """
    Plot the 2pcf multipoles and save the figure if desired.
    
    Parameters
    -----------

    s: numpy.ndarray, required
        The value of s in each bin of the 2pcf.
    
    xi_ell: dict, required
        The 2pcf for each of the desired ells. Individual correlations 
        are stored as arrays.

    filename: str, required
        Location to save the corr file. Please specify
        the full path to the desired fits file.
        
    Returns
    --------

    None

    """
    # Define column lise
    out_cols = []

    # Define s column
    out_cols.append(fits.Column(name='s',array=s,format='E'))
    
    for ell in xi_ell.keys():
        out_cols.append(fits.Column(name=f'xi_l{ell}',
                                    array=xi_ell[ell],format='E'))
        
    # Wrap columns and save file
    out_table = fits.BinTableHDU.from_columns(out_cols)
    out_table.writeto(filename,overwrite=True)

    # No return
    return None

# =========================================================================== #

