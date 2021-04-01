'''
This file gather individual function that are related to external data for the CTA Perseus KSP
'''

import os
import astropy.units as u
import astropy.constants as const
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from minot.model_tools import trapz_loglog


#==================================================
# Get the profile data
#==================================================

def radio_profile_data(distance_correction):
    """
    Get the radio profile data and output a dictionary.
    This is taken from Gitti et al. 2002
    
    Parameters
    ----------
    
    Outputs
    ----------
    - data (dict): a dictionary that contain the data
    
    """
    
    # Extracted from Gitti's paper
    dat_dir    = os.getenv('CTAPHYS_EDAT_DIR')+'Radio'
    prof_file  = dat_dir+'/Perseus_radio_profile_Gitti.txt'
    prof_data11 = pd.read_csv(prof_file, header=None, skiprows=1, index_col=False,
                              names=['radius', 'flux', 'error_m', 'error_p'])
    prof_data1 = {'radius' : prof_data11['radius'].values,
                  'flux'   : 10**prof_data11['flux'],
                  'error_m': np.log(10)*10**prof_data11['flux'].values * prof_data11['error_m'].values / 2,
                  'error_p': np.log(10)*10**prof_data11['flux'].values * prof_data11['error_p'].values / 2}
    
    # Provided by G.F.B.
    prof_data2 = {'radius':(np.array([239.2,217.8,203.5,192.8,185.7,167.8,153.5,
                                      139.2,128.5,110.7,092.8,082.1,067.8,053.6,
                                      042.8,039.3,028.6,014.3]))[::-1], 
                  'flux_m':(np.array([0.003667,0.003667,0.003643,0.007334,0.011000,0.018340,0.036670,
                                      0.036670,0.055010,0.073340,0.073340,0.183400,0.366700,0.550100,
                                      0.733400,0.733400,3.667000,14.66800]))[::-1],
                  'flux_p':(np.array([0.036670,0.036670,0.018340,0.055010,0.055010,0.091680,0.183400,
                                      0.366700,0.366700,0.550100,0.550100,0.733400,0.733400,1.375100,
                                      3.667000,3.667000,7.334100,14.66800]))[::-1]}
    prof_data2['flux']    = 0.5*(prof_data2['flux_m'] + prof_data2['flux_p']) # Mean of flux plus/flux minus
    prof_data2['error_m'] =     (prof_data2['flux']   - prof_data2['flux_m']) # Mean-max difference
    prof_data2['error_p'] =    -(prof_data2['flux']   - prof_data2['flux_p']) # Mean-min difference
    
    # Compute errors
    error_p = np.abs(prof_data2['flux_p']-np.array(prof_data1['flux']))*u.Jy/u.arcmin**2  # 
    error_m = np.abs(prof_data2['flux_m']-np.array(prof_data1['flux']))*u.Jy/u.arcmin**2  # 
    #error   = 0.5 * (prof_data1['error_m']   +   prof_data1['error_p'])*u.Jy/u.arcmin**2  # Symetrized error
    error   = np.array(prof_data1['flux'])*u.Jy/u.arcmin**2 * 0.3                         # Assume 20% error
    
    # Combinaning the two
    prof_data = {'radius':prof_data2['radius']*u.kpc*distance_correction, # standard radius
                 'flux':np.array(prof_data1['flux'])*u.Jy/u.arcmin**2,    # central flux
                 'error':error,                                           # Symmetrixed error
                 'error_p':prof_data1['error_p']*u.Jy/u.arcmin**2,
                 'error_m':prof_data1['error_m']*u.Jy/u.arcmin**2}

    return prof_data


#==================================================
# Get the profile data
#==================================================

def radio_profile_data2(kpcperarcmin):
    """
    Get the radio profile data and output a dictionary
    This is taken from Pedlar et al. 1990
    Assume 10% error

    Parameters
    ----------
    
    Outputs
    ----------
    - data (dict): a dictionary that contain the data
    
    """
    
    # Extracted from Pedlar's paper
    dat_dir    = os.getenv('CTAPHYS_EDAT_DIR')+'Radio'
    prof_file  = dat_dir+'/Perseus_radio_profile_Pedlar1990_1380MHz_beam42x41arcsec.txt'
    prof_data11 = pd.read_csv(prof_file, header=None, skiprows=1, index_col=False,
                              names=['radius', 'flux'])

    sigma = np.sqrt(41*42/60**2)/2.355 # arcmin
    beam = 2*np.pi*sigma**2 # arcmin2/beam

    prof_data1 = {'radius' : prof_data11['radius'].values/60.0*kpcperarcmin,
                  'flux'   : prof_data11['flux'].values/beam,
                  'error'  : prof_data11['flux'].values*0.1/beam,
                  'error_m': prof_data11['flux'].values*0.1/beam,
                  'error_p': prof_data11['flux'].values*0.1/beam}
    
    # Combinaning the two
    prof_data = {'radius':prof_data1['radius']*u.kpc,
                 'flux':prof_data1['flux']*u.Jy/u.arcmin**2,    # central flux
                 'error':prof_data1['error']*u.Jy/u.arcmin**2,   # Symmetrixed error
                 'error_p':prof_data1['error_p']*u.Jy/u.arcmin**2,
                 'error_m':prof_data1['error_m']*u.Jy/u.arcmin**2}

    return prof_data



#==================================================
# Get the spectrum data
#==================================================

def radio_spectrum_data():
    """
    Get the radio spectrum data and output a dictionary
    
    Parameters
    ----------
    
    Outputs
    ----------
    - data (dict): a dictionary that contain the data
    
    """

    spec_data = {'freq':np.array([327,609,1395])*u.MHz, 
                 'flux':np.array([1.2423625254582,0.8798370672098,0.4765784114053]), # This is log F/Jy here
                 'error':np.array([0.12219959266803,0.08961303462322,0.05702647657841])}

    spec_data['error'] = np.log(10)*10**spec_data['flux'] * spec_data['error'] * u.Jy
    spec_data['flux'] = 10**spec_data['flux'] *u.Jy

    return spec_data
    

#==================================================
# Get the spectral profile index data
#==================================================

def radio_index_data(distance_correction):
    """
    Get the radio spectral index profile data and output a 
    dictionary
    
    Parameters
    ----------
    
    Outputs
    ----------
    - data (dict): a dictionary that contain the data
    
    """
    
    idx_data = {'radius': np.array([27,46,70,123,169])*u.kpc*distance_correction, 
                'idx':    np.array([1.10,1.25,1.55,1.9,2.15]),
                'error':  np.array([0.0, 0.125, 0.125, 0.125, 0.125])}

    return idx_data


#==================================================
# Check the radio data consistency
#==================================================

def radio_consistency(radio_data, kpcperarcmin, wspec, check=True):
    """
    Check the radio data consistency between profile
    and spectrum
    
    Parameters
    ----------
    - radio_data (dict): the data
    - kpcperarcmin (float): distance angle conversion
    
    Outputs
    ----------
    flux_correction_factor (float): correction from 
    profile to spectrum
    """

    # Interpolate the data over a requested range
    t_itpl = np.logspace(np.log10((radio_data['info']['spec_Rmin']/kpcperarcmin).to_value('arcmin')),
                         np.log10((radio_data['info']['spec_Rmax']/kpcperarcmin).to_value('arcmin')),
    1000)
    itpl = interp1d((radio_data['profile']['radius']/kpcperarcmin).to_value('arcmin'), 
                    radio_data['profile']['flux'].to_value('Jy arcmin-2'), 
                    kind='linear', fill_value='extrapolate')
    p_itpl = itpl(t_itpl)

    # Compute the flux expected from the profile
    flux_from_profile = trapz_loglog(2*np.pi*t_itpl*p_itpl, t_itpl)*u.Jy
    flux_correction_factor = flux_from_profile/radio_data['spectrum']['flux'][wspec]
    
    # Check the interpolation
    if check:
        import matplotlib.pyplot as plt
        plt.loglog((radio_data['profile']['radius']/kpcperarcmin).to_value('arcmin'), 
                   radio_data['profile']['flux'].to_value('Jy arcmin-2'))
        plt.loglog(t_itpl, p_itpl)
        plt.xlabel('radius (arcmin)')
        plt.ylabel('Surface brightness (Jy arcmin-2)')
        print('Flux from profile:', flux_from_profile)
        print('Flux from spectrum:', radio_data['info']['prof_freq'],
              radio_data['spectrum']['freq'][wspec], radio_data['spectrum']['flux'][wspec])
        print('flux correcction (to be applied to profile)', flux_correction_factor)

    return flux_correction_factor


#==================================================
# Get the radio data
#==================================================

def get_radio_data(cosmo, redshift, prof_file='Gitti2002'):
    """
    Get the total radio data in one single dictionary
    
    Parameters
    ----------
    - cosmo (astropy.cosmology): a cosmology object
    
    Outputs
    ----------
    radio_data (dict): dictionary containing the radio data
    """

    # Cosmology
    distance_correction = 50.0/cosmo.H0.to_value('km s-1 Mpc-1')
    kpcperarcmin = cosmo.kpc_proper_per_arcmin(redshift)
    
    # Information about the data
    if prof_file == 'Gitti2002':
        info = {'spec_Rmin' : 1*u.kpc,                       # Radius down to which the flux is integrated
                'spec_Rmax' : 15.0/2*u.arcmin*kpcperarcmin,  # Radius up to which the flux is integrated
                'prof_Rmin' : 30*u.kpc*distance_correction,  # Radius down to which the model ok (due to NGC1275)
                'prof_Rmax' : 500*u.kpc*distance_correction, # Radius down to which the model ok (due to NGC1275)
                'prof_freq' : 327*u.MHz,                     # Frequency at which the profile is extracted
                'idx_freq1' : 327*u.MHz,                     # Start frequency for spectral index calculation
                'idx_freq2' : 609*u.MHz,                     # End frequency for spectral index calculation      
                'idx_Rmin'  : 30*u.kpc*distance_correction,  # Radius down to which the model ok (due to NGC1275)
                'idx_Rmax'  : 500*u.kpc*distance_correction} # Radius down to which the model ok (due to NGC1275)
    if prof_file == 'Pedlar1990':
        info = {'spec_Rmin' : 1.0*u.kpc,                     # Radius down to which the flux is integrated
                'spec_Rmax' : 15.0/2*u.arcmin*kpcperarcmin,  # Radius up to which the flux is integrated
                'prof_Rmin' : 23*u.kpc,                      # Radius down to which the model ok (due to NGC1275)
                'prof_Rmax' : 80*u.kpc,                      # Radius down to which the model ok (due to NGC1275)
                'prof_freq' : 1380*u.MHz,                    # Frequency at which the profile is extracted
                'idx_freq1' : 327*u.MHz,                     # Start frequency for spectral index calculation
                'idx_freq2' : 609*u.MHz,                     # End frequency for spectral index calculation      
                'idx_Rmin'  : 30*u.kpc*distance_correction,  # Radius down to which the model ok (due to NGC1275)
                'idx_Rmax'  : 500*u.kpc*distance_correction} # Radius down to which the model ok (due to NGC1275)

    # Data
    if prof_file == 'Gitti2002':
        prof_data = radio_profile_data(distance_correction)
    if prof_file == 'Pedlar1990':        
        prof_data = radio_profile_data2(kpcperarcmin.value)

    spec_data = radio_spectrum_data()
    idx_data  = radio_index_data(distance_correction)
    
    radio_data = {'info':     info,
                  'profile' : prof_data,
                  'spectrum': spec_data,
                  'index':    idx_data}

    # Profile - spectrum cross calibration
    if prof_file == 'Gitti2002': wspec = 0
    if prof_file == 'Pedlar1990': wspec = 2
    flux_correction = radio_consistency(radio_data, kpcperarcmin, wspec=wspec, check=True)
    print(flux_correction)
    if prof_file == 'Gitti2002':
        print('---> Apply flux correction for Gitti 2002 data')
        radio_data['profile']['flux'] = radio_data['profile']['flux']/flux_correction
        radio_data['profile']['error'] = radio_data['profile']['error']/flux_correction
        radio_data['profile']['error_p'] = radio_data['profile']['error_p']/flux_correction
        radio_data['profile']['error_m'] = radio_data['profile']['error_m']/flux_correction
    if prof_file == 'Pedlar1990': 
        print('---> Check spectrum/profile consistency for Pedlar 1990 data, but no correction.')

    return radio_data
