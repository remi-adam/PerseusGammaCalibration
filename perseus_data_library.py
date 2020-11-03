'''
This file gather individual function that are related to external data for the CTA Perseus KSP
'''

import os
import astropy.units as u
import astropy.constants as const

dat_dir = os.getenv('CTAPHYS_EDAT_DIR')+'Radio'


#==================================================
# Get the profile data
#==================================================

def radio_profile_data():
    """
    Get the radio profile data and output a dictionary
    
    Parameters
    ----------
    
    Outputs
    ----------
    - data (dict): a dictionary that contain the data
    
    """

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

#==================================================
# Get the spectral profile index data
#==================================================

def radio_index_data():
    """
    Get the radio spectral index profile data and output a 
    dictionary
    
    Parameters
    ----------
    
    Outputs
    ----------
    - data (dict): a dictionary that contain the data
    
    """



#==================================================
# Check the radio data consistency
#==================================================

def radio_consistency():
    """
    Check the radio data consistency between profile
    and spectrum
    
    Parameters
    ----------
    
    Outputs
    ----------
    flux_correction_factor (float): correction from 
    profile to spectrum
    """
