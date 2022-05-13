'''
This file gather individual function that define the Perseus cluster model for the CTA KSP
'''

import numpy as np
import os
import copy
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(70, 0.3, m_nu=[0.0,0.0,0.06]*u.eV, Tcmb0=2.725, Ob0=0.0486)
import minot


#==================================================
# Define the default cluster model
#==================================================

def default_model(directory=None):
    """
    Define the default cluster model
    
    Parameters
    ----------
    
    Outputs
    ----------
    - cluster (minot object): 
    
    """

    if directory is None:
        outdir = os.getenv('CTAPHYS_OUT_DIR')+'Perseus_KSP_calibration'
    else:
        outdir = directory
        
    #---------- Global parameters
    NR500_trunc = 3
    redshift    = 0.017284         # Hitomi 2018
    theta500    = 59.7*u.arcmin    # Urban (2014) PIP V UPP best fit
    RA          = 49.950667*u.deg  # NED (NGC 1275)
    Dec         = 41.511696*u.deg  # NED (NGC 1275)

    Yhe         = 0.2735           # Baseline in minot
    Abundance   = 0.3              # In agreement with Werner et al. (2013)

    EBL         = 'dominguez'      # Arbitrary

    #Epmin = (1+0.8**2)**0.5*(const.m_p*const.c**2)).to('GeV') # Aleksic et al. (2012)
    Epmin = minot.ClusterTools.cluster_spectra.pp_pion_kinematic_energy_threshold()*u.GeV
    Epmax = 10*u.PeV
    
    #---------- Define cluster
    cluster = minot.Cluster(name='Perseus',
                            redshift=redshift,
                            RA=RA, Dec=Dec,
                            cosmology=cosmo,
                            silent=True,
                            output_dir=outdir)
    cluster.theta500     = theta500
    cluster.R_truncation = NR500_trunc*cluster.R500

    cluster.helium_mass_fraction = Yhe
    cluster.abundance            = Abundance

    cluster.EBL_model = EBL
    cluster.Epmin     = Epmin
    cluster.Epmax     = Epmax

    #---------- Set the model
    cluster = set_thermal_model(cluster)
    cluster = set_magnetic_field_model(cluster, case='Taylor2006')
    cluster = set_pure_hadronic_model(cluster, ('density', 1.0), 1e-2, 2.3)

    return cluster


#==================================================
# Set the gas thermal model to the cluster
#==================================================

def set_thermal_model(cluster_in):
    """
    Set the thermal gas model
    
    Parameters
    ----------
    - cluster_in (minot object): input cluster model
    
    Outputs
    ----------
    - cluster_out (minot object): modified cluster model
    
    """

    cluster_out = copy.deepcopy(cluster_in)

    # Density
    cor1 = cluster_out.cosmo.H0.value / 50.0 # From Churazov 2003
    cor2 = cluster_out.cosmo.H0.value / 70.0 # From Urban 2014
    arcmin2kpc = cluster_out.cosmo.kpc_proper_per_arcmin(cluster_out.redshift)
    
    cluster_out.density_gas_model = {'name':'doublebeta',
                                     'beta1':1.2,
                                     'r_c1':80*u.kpc * cor1**-1,
                                     'n_01':0.039*u.cm**-3 * cor1**0.5,
                                     'beta2':0.71,
                                     'r_c2':(13.18*u.arcmin*arcmin2kpc).to('kpc'),
                                     'n_02':0.0036*u.cm**-3 * cor2**0.5}
    
    # Pressure
    radius = np.logspace(-1,5,10000)
    cor1 = cluster_out.cosmo.H0.value / 50.0 # From Churazov 2003
    cor2 = cluster_out.cosmo.H0.value / 70.0 # From Urban 2014
    T_e = 7*u.keV*(1+(radius/100*cor1)**3)/(2.3+(radius/100*cor1)**3)*(1+(radius/1600*cor2)**1.7)**-1
    n_e = cluster_out.get_density_gas_profile(radius*u.kpc)[1]
    P_e = n_e*T_e
    cluster_out.pressure_gas_model = {'name':'User', 'radius':radius*u.kpc, 'profile':P_e}
    
    return cluster_out


#==================================================
# Set the magnetic field model to the cluster
#==================================================

def set_magnetic_field_model(cluster_in, case='Taylor2006'):
    """
    Set the magnetic field model
    
    Parameters
    ----------
    - cluster_in (minot object): input cluster model
    
    Outputs
    ----------
    - cluster_out (minot object): modified cluster model
    - case (string): which model to use
    
    """

    cluster_out = copy.deepcopy(cluster_in)

    #----- Useful general quantities
    radius = np.logspace(0,4,1000)*u.kpc
    n0_coma = 2.89e-3*u.cm**-3 * (50.0/cluster_out.cosmo.H0.value)**-0.5

    #----- Cases
    if case == 'Taylor2006':
        cluster_out.name = r'Taylor (2006) + $\eta=2/3$'
        cluster_out.set_magfield_isodens_scal_param(Bnorm=1*u.uG, r0=10*u.kpc, scal=1.0)
        r, B = cluster_out.get_magfield_profile(radius)
        r_ref, B_ref = cluster_out.get_magfield_profile(10*u.kpc)
        cluster_out.magfield_model = {'name':'User', 'radius':radius, 'profile':25*u.uG*(B/B_ref[0])**(2.0/3)}

    elif case == 'Walker2017':
        cluster_out.name = 'Walker (2017) + pressure model'
        r, P = cluster_out.get_pressure_gas_profile(radius)
        Y = cluster_out.helium_mass_fraction
        Z = cluster_out.metallicity_sol*cluster_out.abundance
        mu_gas,mu_e,mu_p,mu_alpha = minot.ClusterTools.cluster_global.mean_molecular_weight(Y=Y, Z=Z)
        B = (2*const.mu0*(mu_e/mu_gas)*P/200)**0.5
        cluster_out.magfield_model = {'name':'User', 'radius':radius, 'profile':B}
        
    elif case == 'Bonafede2010best':
        cluster_out.name = r'Bonafede (2010, Coma) + best $\eta$ scaling'
        n_e = cluster_out.get_density_gas_profile(radius)[1]
        cluster_out.magfield_model = {'name':'User', 'radius':radius,
                                      'profile':4.7*u.uG * ((n_e/n0_coma).to_value(''))**0.5}
        
    elif case == 'Bonafede2010low':
        cluster_out.name = r'Bonafede (2010, Coma) + lower $\eta$ scaling'
        n_e = cluster_out.get_density_gas_profile(radius)[1]
        cluster_out.magfield_model = {'name':'User', 'radius':radius,
                                      'profile':3.9*u.uG * ((n_e/n0_coma).to_value(''))**0.4}
    elif case == 'Bonafede2010up':
        cluster_out.name = r'Bonafede (2010, Coma) + upper $\eta$ scaling'
        n_e = cluster_out.get_density_gas_profile(radius)[1]
        cluster_out.magfield_model = {'name':'User', 'radius':radius,
                                      'profile':5.4*u.uG * ((n_e/n0_coma).to_value(''))**0.9}
    elif case == 'Bonafede2010std':
        cluster_out.name = r'Bonafede (2010, Coma) + $\eta=2/3$ scaling'
        n_e = cluster_out.get_density_gas_profile(radius)[1]
        cluster_out.magfield_model = {'name':'User', 'radius':radius,
                                      'profile':5.0*u.uG * ((n_e/n0_coma).to_value(''))**(2.0/3)}
    else:
        print('Available models are: ')
        print('   Taylor2006, Walker2017, Bonafede2010best, Bonafede2010low, Bonafede2010up, Bonafede2010std')
        print('Doing nothing')

    return cluster_out


#==================================================
# Set the pure hadronic CR model to the cluster
#==================================================

def set_pure_hadronic_model(cluster_in, scaling, Xcrp, slope):
    """
    Set the CR to pure hadronic model
    
    Parameters
    ----------
    - cluster_in (minot object): input cluster model
    - scaling (tupple): scaling[0] is the reference thermodynamic 
    quantity ('density' or 'pressure') and scaling[1] is the scaling 
    value eta
    - Xcrp (float): amplitude as U_CRp / U_th at R500
    - slope (float): the slope for the CRp

    Outputs
    ----------
    - cluster_out (minot object): modified cluster model
    
    """
    
    cluster_out = copy.deepcopy(cluster_in)

    cluster_out.X_cre1_E = {'X':0.00, 'R_norm':cluster_out.R500}
    
    cluster_out.X_crp_E  = {'X':Xcrp, 'R_norm':cluster_out.R500}
    #cluster_out.spectrum_crp_model = {'name':'PowerLaw', 'Index':slope}
    cluster_out.spectrum_crp_model = {'name':'MomentumPowerLaw', 'Index':slope, 'Mass':(const.m_p*const.c**2)}
    radius = np.logspace(0,4,1000)*u.kpc

    if scaling[0] == 'density':
        if cluster_out.density_gas_model['name'] in ['doublebeta', 'User']:
            n_e = cluster_out.get_density_gas_profile(radius)[1]
            cluster_out.density_crp_model = {'name':'User', 'radius':radius,
                                             'profile':n_e.value**scaling[1]}
            
        else:
            cluster_out.set_density_crp_isodens_scal_param(scal=scaling[1])
        
    elif scaling[0] == 'pressure':
        if cluster_out.pressure_gas_model['name'] in ['doublebeta', 'User']:
            p_e = cluster_out.get_pressure_gas_profile(radius)[1]
            cluster_out.density_crp_model = {'name':'User', 'radius':radius,
                                             'profile':p_e.value**scaling[1]}

        else:        
            cluster_out.set_density_crp_isobaric_scal_param(scal=scaling[1])
        
    else:
        raise ValueError('scaling[0] should be "density" or "pressure"')
        
    return cluster_out


#==================================================
# Set the pure leptonic CR model to the cluster
#==================================================

def set_pure_leptonic_model(cluster_in, scaling, Xcre1, slope):
    """
    Set the CR to pure leptonic model: power law injection + steady losses
    
    Parameters
    ----------
    - cluster_in (minot object): input cluster model
    - scaling (tupple): scaling[0] is the reference thermodynamic 
    quantity ('density' or 'pressure') and scaling[1] is the scaling 
    value eta
    - Xcre (float): amplitude as U_CRe / U_th at R500
    - slope (float): the slope for the CRp
    - Ecut (quantity): energy cutoff, homogeneous to GeV

    Outputs
    ----------
    - cluster_out (minot object): modified cluster model
    
    """
    
    cluster_out = copy.deepcopy(cluster_in)

    cluster_out.X_crp_E  = {'X':0.00, 'R_norm':cluster_out.R500}
    
    cluster_out.X_cre1_E = {'X':Xcre1, 'R_norm':cluster_out.R500}
    cluster_out.spectrum_cre1_model = {'name': 'PowerLaw',
                                       'Index': slope}
    cluster_out.cre1_loss_model = 'Steady'
    
    radius = np.logspace(0,4,1000)*u.kpc

    if scaling[0] == 'density':
        if cluster_out.density_gas_model['name'] in ['doublebeta', 'User']:
            n_e = cluster_out.get_density_gas_profile(radius)[1]
            cluster_out.density_cre1_model = {'name':'User', 'radius':radius,
                                              'profile':n_e.value**scaling[1]}
            
        else:
            cluster_out.set_density_cre1_isodens_scal_param(scal=scaling[1])
        
    elif scaling[0] == 'pressure':
        if cluster_out.pressure_gas_model['name'] in ['doublebeta', 'User']:
            p_e = cluster_out.get_pressure_gas_profile(radius)[1]
            cluster_out.density_cre1_model = {'name':'User', 'radius':radius,
                                              'profile':p_e.value**scaling[1]}

        else:        
            cluster_out.set_density_cre1_isobaric_scal_param(scal=scaling[1])
        
    else:
        raise ValueError('scaling[0] should be "density" or "pressure"')
        
    return cluster_out
