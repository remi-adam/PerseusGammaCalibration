"""
This script is part of the CTA key science project dedicated to Perseus. 
Here, we constrain the Perseus CR distribution (spectral and spatial) using 
radio data, via an MCMC approach. This includes both pure hadronic and pure 
leptonioc models.
"""

#========================================
# Imports
#========================================

import matplotlib
matplotlib.use('Agg')

from multiprocessing import Pool, cpu_count
import warnings
import astropy.units as u
import numpy as np
import pandas as pd
import emcee
import corner
import matplotlib.pyplot as plt
import copy
import pickle
import os
from scipy.optimize import curve_fit

from PerseusGammaCalibration import perseus_model_library
from PerseusGammaCalibration import perseus_data_library
from kesacco.Tools.plotting import seaborn_corner
from kesacco.Tools import mcmc_common
from kesacco.Tools import plotting
import minot

dict_base = {'font.size':        16,
             'legend.fontsize':  16,
             'xtick.labelsize':  16,
             'ytick.labelsize':  16,
             'axes.labelsize':   16,
             'axes.titlesize':   16,
             'figure.titlesize': 16,    
             'figure.figsize':[8.0, 6.0],
             'figure.subplot.right':0.97,
             'figure.subplot.left':0.15,
             'font.family':'serif',
             'figure.facecolor': 'white',
             'legend.frameon': True}
plt.rcParams.update(dict_base)

#========================================
# MCMC post analysys
#========================================

def post_analysis(cluster, radio_data, param_name, par_min, par_max, burnin,
                  conf=68.0, Nmc=100, model_case='Hadronic'):
    '''
    Extract statistics and plots after runing the MCMC

    Parameters
    ----------
    - cluster (minot object): cluster object
    - radio_data (dict): radio data
    - param_name (list): list of the names of the parameters
    - par_min/max (list): list of min and max value for the parameters
    - burnin (int): the burn-in phase of the chain
    - conf (float): the confidence limit interval
    - Nmc (int): the number of Monte Carlo sample to draw

    Output
    ------
    - Plots and statistics
    '''

    #========== Get the chains    
    #----- Restore chains
    with open(cluster.output_dir+'/'+model_case+'_sampler.pkl', 'rb') as f:
        sampler = pickle.load(f)

    #----- Plot acceptance
    # Acceptance
    sample_acc = (np.roll(sampler.lnprobability,1,axis=1) != sampler.lnprobability)
    sample_acc = sample_acc[:,1:]

    acc_mean   = np.sum(sample_acc)/float(len(sample_acc.flatten()))
    acc_sample = np.sum(sample_acc, axis=0)/float(len(sample_acc[:,0]))
    acc_chain  = np.sum(sample_acc, axis=1)/float(len(sample_acc[0,:]))
    
    plt.figure(10)
    plt.plot(sampler.acceptance_fraction, 'k')
    plt.plot(acc_chain, 'r', linestyle='--')
    plt.ylim(0,1)
    plt.xlim(0,len(sampler.acceptance_fraction)-1)
    plt.xlabel('Chain number')
    plt.ylabel('Mean sample acceptance')
    plt.savefig(cluster.output_dir+'/'+model_case+'_chains_acceptance_chain.pdf')
    plt.close()
    
    plt.figure(11)
    plt.plot(acc_sample, 'k')
    plt.ylim(0,1)
    plt.xlim(0,len(acc_sample)-1)
    plt.xlabel('Sample number')
    plt.ylabel('Mean chain acceptance')
    plt.savefig(cluster.output_dir+'/'+model_case+'_chains_acceptance_sample.pdf')
    plt.close()
    print('  --> Acceptence plot done')
    
    #----- Burn in
    param_chains = sampler.chain[:, burnin:, :]
    lnL_chains = sampler.lnprobability[:, burnin:]
    ndim = param_chains.shape[2]
    print('  --> Burn-in done')

    #----- Chain covariance
    chain_list = []
    for i in range(ndim):
        chain_list.append(param_chains[:,:,i].flatten())
    par_cov = np.cov(np.array(chain_list))
    print('---- Param chain covariance^1/2:')
    print(par_cov**0.5)
    print('')
    print('  --> Chain covariance done')
        
    #----- Get the best fit parameters
    wbest = (lnL_chains == np.amax(lnL_chains))
    param_best = []
    for i in range(ndim):
        param_best.append(((param_chains[:,:,i])[wbest])[0])
    print('  --> Best-fit parameters done')

    #----- MC parameters
    param_flat = param_chains.reshape(param_chains.shape[0]*param_chains.shape[1],param_chains.shape[2])
    Nsample = len(param_flat[:,0])-1
    param_MC = np.zeros((Nmc, ndim))
    for i in range(Nmc):
        param_MC[i,:] = param_flat[np.random.randint(0, high=Nsample), :] # randomly taken from chains
    print('  --> MC parameters done')
        
    #========== Statistics results
    #----- Chain statistics
    chainfname = cluster.output_dir+'/'+model_case+'_chainstat.txt'
    par_best, par_percentile = mcmc_common.chains_statistics(param_chains, lnL_chains,
                                                             parname=param_name, conf=conf, show=True,
                                                             outfile=chainfname)
    print('  --> Chain statistics done')
        
    #----- Parameter space
    try:
        mcmc_common.chains_plots(param_chains, param_name, cluster.output_dir+'/'+model_case+'_chains',
                                 par_best=par_best, par_percentile=par_percentile, conf=conf,
                                 par_min=par_min, par_max=par_max)
    except:
        print('Error in chains_plots --> skip it')
    print('  --> Chain plots done')

    #========== Data versus model
    # Best-fit
    if model_case == 'Hadronic':
        prof_best, spec_best, idx_best = model_hadronic(par_best, cluster, radio_data)
    if model_case == 'Leptonic':
        prof_best, spec_best, idx_best = model_leptonic(par_best, cluster, radio_data)
                
    # MC sampling
    prof_mc = []
    spec_mc = []
    idx_mc  = []
    
    for imc in range(Nmc):
        if model_case == 'Hadronic':
            prof_mci, spec_mci, idx_mci = model_hadronic(param_MC[imc,:], cluster, radio_data)
        if model_case == 'Leptonic':
            prof_mci, spec_mci, idx_mci = model_leptonic(param_MC[imc,:], cluster, radio_data)            
        prof_mc.append(prof_mci)
        spec_mc.append(spec_mci)
        idx_mc.append(idx_mci)

    # Limits
    prof_u = np.percentile(np.array(prof_mc), 100-(100-conf)/2.0, axis=0)*prof_mc[0].unit
    prof_d = np.percentile(np.array(prof_mc), (100-conf)/2.0, axis=0)*prof_mc[0].unit
    
    spec_u = np.percentile(np.array(spec_mc), 100-(100-conf)/2.0, axis=0)*spec_mc[0].unit
    spec_d = np.percentile(np.array(spec_mc), (100-conf)/2.0, axis=0)*spec_mc[0].unit
    
    idx_u = np.percentile(np.array(idx_mc), 100-(100-conf)/2.0, axis=0)
    idx_d = np.percentile(np.array(idx_mc), (100-conf)/2.0, axis=0)

    #----- Spectrum
    fig = plt.figure(0, figsize=(8, 6))
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_mc[imc].to_value('Jy'),
                 color='blue', alpha=0.05)
        
    # Limits
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_u.to_value('Jy'),
             color='blue', linewidth=2, linestyle='--', label=str(conf)+'% C.L.')
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_d.to_value('Jy'),
             color='blue', linewidth=2, linestyle='--')
    plt.fill_between(radio_data['spectrum']['freq'].to_value('MHz'), spec_d.to_value('Jy'),
                     spec_u.to_value('Jy'), color='blue', alpha=0.2)
    
    # Best fit and data
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_best.to_value('Jy'),
             color='blue', linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['spectrum']['freq'].to_value('MHz'), radio_data['spectrum']['flux'].to_value('Jy'),
                 radio_data['spectrum']['error'].to_value('Jy'),
                 marker='o', color='k', linestyle='', label='Data')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('frequency (MHz)')
    plt.ylabel('flux (Jy)')
    plt.legend()
    plt.savefig(cluster.output_dir+'/'+model_case+'_Radio_spectrum.pdf')
    plt.close()

    #----- Profile
    fig = plt.figure(0, figsize=(8, 6))
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_mc[imc].to_value('Jy arcmin-2'),
                 color='blue', alpha=0.05)
        
    # Limits    
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_u.to_value('Jy arcmin-2'),
             color='blue', linewidth=2, linestyle='--', label=str(conf)+'% C.L.')
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_d.to_value('Jy arcmin-2'),
             color='blue', linewidth=2, linestyle='--')
    plt.fill_between(radio_data['profile']['radius'].to_value('kpc'), prof_d.to_value('Jy arcmin-2'),
                     prof_u.to_value('Jy arcmin-2'), color='blue', alpha=0.2)
    
    # Best and data
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_best.to_value('Jy arcmin-2'),
             color='blue', linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['profile']['radius'].to_value('kpc'),
                 radio_data['profile']['flux'].to_value('Jy arcmin-2'), 
                 yerr=radio_data['profile']['error'].to_value('Jy arcmin-2'),
                 marker='o', linestyle='', color='k', label='Data')
    
    plt.plot([radio_data['info']['prof_Rmin'].to_value('kpc'),radio_data['info']['prof_Rmin'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.plot([radio_data['info']['prof_Rmax'].to_value('kpc'),radio_data['info']['prof_Rmax'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.fill_between([0,radio_data['info']['prof_Rmin'].to_value('kpc')], [0,0], [1e6,1e6],
                     color='grey', alpha=0.2, label='Excluded region')
    plt.fill_between([radio_data['info']['prof_Rmax'].to_value('kpc'),
                      radio_data['info']['prof_Rmax'].to_value('kpc')*1e6], [0,0], [1e6,1e6],
                     color='grey', alpha=0.2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('radius (kpc)')
    plt.ylabel('surface brightness (Jy/arcmin$^2$)')
    if basedata == 'Gitti2002':
        plt.xlim(10,250)
        plt.ylim(1e-3,1e1)
    if basedata == 'Pedlar1990':
        plt.xlim(13,100)
        plt.ylim(5e-3,5e-1)
    plt.savefig(cluster.output_dir+'/'+model_case+'_Radio_profile.pdf')
    plt.close()
    
    #----- Spectral index
    if fit_index:
        fig = plt.figure(0, figsize=(8, 6))
        # MC
        for imc in range(Nmc):
            plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_mc[imc], color='blue', alpha=0.05)
            
        # Limits
        plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_u, color='blue',
                 linewidth=2, linestyle='--', label=str(conf)+'% C.L.')
        plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_d, color='blue',
                 linewidth=2, linestyle='--')
        plt.fill_between(radio_data['index']['radius'].to_value('kpc'), idx_d, idx_u,
                         color='blue', alpha=0.2)
        
        # Best and data    
        plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_best, color='blue',
                 linewidth=3, label='Best-fit')
        plt.errorbar(radio_data['index']['radius'].to_value('kpc'), radio_data['index']['idx'],
                     yerr=radio_data['index']['error'],
                     marker='o', linestyle='', color='k', label='Data')
        
        plt.plot([radio_data['info']['idx_Rmin'].to_value('kpc'),radio_data['info']['idx_Rmin'].to_value('kpc')],
                 [0,1e6], linestyle='--', color='grey')
        plt.plot([radio_data['info']['idx_Rmax'].to_value('kpc'),radio_data['info']['idx_Rmax'].to_value('kpc')],
                 [0,1e6], linestyle='--', color='grey')
        plt.fill_between([0,radio_data['info']['idx_Rmin'].to_value('kpc')], [0,0], [1e6,1e6], color='grey',
                         alpha=0.1, label='Excluded region')
        plt.fill_between([radio_data['info']['idx_Rmax'].to_value('kpc'), np.inf], [0,0], [1e6,1e6],
                         color='grey', alpha=0.1)
        plt.xscale('log')
        plt.yscale('linear')
        plt.xlabel('radius (kpc)')
        plt.ylabel('spectral index')
        plt.xlim(10,250)
        plt.ylim(0.5,2.5)
        plt.savefig(cluster.output_dir+'/'+model_case+'_Radio_index.pdf')
        plt.close()

    print('  --> Radio plots done')
    
    #========== Implication for gamma rays
    energy = np.logspace(-2,6,100)*u.GeV
    radius = np.logspace(0,4,100)*u.kpc
    
    #---------- Best-fit
    if model_case == 'Hadronic':
        cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', param_best[1]),
                                                                param_best[0]*1e-2, param_best[2])
    if model_case == 'Leptonic':
        cluster = perseus_model_library.set_pure_leptonic_model(cluster, ('density', param_best[1]),
                                                                param_best[0]*1e-5, param_best[2])
        if app_steady: cluster.cre1_loss_model = 'Steady'

    # Hadronic
    E, dN_dEdSdt = cluster.get_gamma_spectrum(energy, 
                                              Rmin=None, Rmax=cluster.R500,
                                              type_integral='cylindrical',
                                              Rmin_los=None, NR500_los=5.0)
    
    r, dN_dSdtdO = cluster.get_gamma_profile(radius, 
                                             Emin=150*u.GeV, Emax=50*u.TeV, 
                                             Energy_density=False, Rmin_los=None, NR500_los=5.0)

    fluxH = cluster.get_gamma_flux(Emin=150*u.GeV, Emax=None, Energy_density=False,
                                   Rmax=cluster.R500, type_integral='cylindrical').to_value('s-1 cm-2')
    
    # Inverse Compton
    E, dNIC_dEdSdt = cluster.get_ic_spectrum(energy, 
                                             Rmin=None, Rmax=cluster.R500,
                                             type_integral='cylindrical',
                                             Rmin_los=None, NR500_los=5.0)
    
    r, dNIC_dSdtdO = cluster.get_ic_profile(radius, 
                                            Emin=150*u.GeV, Emax=50*u.TeV, 
                                            Energy_density=False, Rmin_los=None, NR500_los=5.0)

    fluxIC = cluster.get_ic_flux(Emin=150*u.GeV, Emax=None, Energy_density=False,
                                 Rmax=cluster.R500, type_integral='cylindrical').to_value('s-1 cm-2')

    #---------- Monte Carlo sampling: hadronic
    prof_g_mc = []
    spec_g_mc = []
    flux_g_mc = []
    
    for imc in range(Nmc):
        if model_case == 'Hadronic':
            cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', param_MC[imc,1]), 
                                                                    param_MC[imc,0]*1e-2, param_MC[imc,2])
        if model_case == 'Leptonic':
            cluster = perseus_model_library.set_pure_leptonic_model(cluster, ('density', param_MC[imc,1]),
                                                                    param_MC[imc,0]*1e-5, param_MC[imc,2])
            if app_steady: cluster.cre1_loss_model = 'Steady'

        spec_g_mci = cluster.get_gamma_spectrum(energy, Rmin=None, Rmax=cluster.R500,
                                                type_integral='cylindrical',
                                                Rmin_los=None, NR500_los=5.0)[1]
    
        prof_g_mci = cluster.get_gamma_profile(radius, Emin=150*u.GeV, Emax=50*u.TeV, 
                                               Energy_density=False, Rmin_los=None, NR500_los=5.0)[1]

        flux_g_mci = cluster.get_gamma_flux(Emin=150*u.GeV, Emax=None, Energy_density=False,
                                            Rmax=cluster.R500, type_integral='cylindrical').to_value('s-1 cm-2')
    
        prof_g_mc.append(prof_g_mci)
        spec_g_mc.append(spec_g_mci)
        flux_g_mc.append(flux_g_mci)
    
    #---------- Limits: hadronic
    dN_dEdSdt_u = np.percentile(np.array(spec_g_mc), 100-(100-conf)/2.0, axis=0)*spec_g_mc[0].unit
    dN_dEdSdt_d = np.percentile(np.array(spec_g_mc), (100-conf)/2.0, axis=0)*spec_g_mc[0].unit
    
    dN_dSdtdO_u = np.percentile(np.array(prof_g_mc), 100-(100-conf)/2.0, axis=0)*prof_g_mc[0].unit
    dN_dSdtdO_d = np.percentile(np.array(prof_g_mc), (100-conf)/2.0, axis=0)*prof_g_mc[0].unit

    fluxH_u = np.percentile(np.array(flux_g_mc), 100-(100-conf)/2.0)
    fluxH_d = np.percentile(np.array(flux_g_mc), (100-conf)/2.0)
    
    #---------- Monte Carlo sampling: IC
    prof_ic_mc = []
    spec_ic_mc = []
    flux_ic_mc = []
    
    for imc in range(Nmc):
        if model_case == 'Hadronic':
            cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', param_MC[imc,1]), 
                                                                    param_MC[imc,0]*1e-2, param_MC[imc,2])
        if model_case == 'Leptonic':
            cluster = perseus_model_library.set_pure_leptonic_model(cluster, ('density', param_MC[imc,1]),
                                                                    param_MC[imc,0]*1e-5, param_MC[imc,2])
        if app_steady: cluster.cre1_loss_model = 'Steady'

        spec_ic_mci = cluster.get_ic_spectrum(energy, Rmin=None, Rmax=cluster.R500,
                                              type_integral='cylindrical',
                                              Rmin_los=None, NR500_los=5.0)[1]
    
        prof_ic_mci = cluster.get_ic_profile(radius, Emin=150*u.GeV, Emax=50*u.TeV, 
                                             Energy_density=False, Rmin_los=None, NR500_los=5.0)[1]

        flux_ic_mci = cluster.get_ic_flux(Emin=150*u.GeV, Emax=None, Energy_density=False,
                                          Rmax=cluster.R500, type_integral='cylindrical').to_value('s-1 cm-2')

        prof_ic_mc.append(prof_ic_mci)
        spec_ic_mc.append(spec_ic_mci)
        flux_ic_mc.append(flux_ic_mci)
    
    #---------- Limits: IC
    dNIC_dEdSdt_u = np.percentile(np.array(spec_ic_mc), 100-(100-conf)/2.0, axis=0)*spec_ic_mc[0].unit
    dNIC_dEdSdt_d = np.percentile(np.array(spec_ic_mc), (100-conf)/2.0, axis=0)*spec_ic_mc[0].unit
    
    dNIC_dSdtdO_u = np.percentile(np.array(prof_ic_mc), 100-(100-conf)/2.0, axis=0)*prof_ic_mc[0].unit
    dNIC_dSdtdO_d = np.percentile(np.array(prof_ic_mc), (100-conf)/2.0, axis=0)*prof_ic_mc[0].unit

    fluxIC_u = np.percentile(np.array(flux_ic_mc), 100-(100-conf)/2.0)
    fluxIC_d = np.percentile(np.array(flux_ic_mc), (100-conf)/2.0)
    
    #========== Figure
    base_cluster = perseus_model_library.default_model(directory=cluster.output_dir)
    base_cluster = perseus_model_library.set_pure_hadronic_model(base_cluster, ('density', 1.0), 1e-2, 2.3)
    
    bid, base_dN_dEdSdt = base_cluster.get_gamma_spectrum(energy, 
                                                          Rmin=None, Rmax=base_cluster.R500,
                                                          type_integral='cylindrical',
                                                          Rmin_los=None, NR500_los=5.0)
    bid, base_dN_dSdtdO = base_cluster.get_gamma_profile(radius, 
                                                         Emin=150*u.GeV, Emax=50*u.TeV, 
                                                         Energy_density=False, Rmin_los=None, NR500_los=5.0)
    
    bid, base_dNIC_dEdSdt = base_cluster.get_ic_spectrum(energy, 
                                                         Rmin=None, Rmax=base_cluster.R500,
                                                         type_integral='cylindrical',
                                                         Rmin_los=None, NR500_los=5.0)
    
    bid, base_dNIC_dSdtdO = base_cluster.get_ic_profile(radius, 
                                                        Emin=150*u.GeV, Emax=50*u.TeV, 
                                                        Energy_density=False, Rmin_los=None, NR500_los=5.0)
    
    #----- Spectrum
    fig = plt.figure(0, figsize=(8, 6))
    # MC
    for imc in range(Nmc):
        if imc == 0:
            plt.plot(E.to_value('GeV'), (E**2*spec_g_mc[imc]).to_value('MeV cm-2 s-1'), color='blue',
                     alpha=0.05, label='Monte Carlo')
        else:
            plt.plot(E.to_value('GeV'), (E**2*spec_g_mc[imc]).to_value('MeV cm-2 s-1'),
                     color='blue', alpha=0.05)
    for imc in range(Nmc):
        plt.plot(E.to_value('GeV'), (E**2*spec_ic_mc[imc]).to_value('MeV cm-2 s-1'),
                 color='magenta', alpha=0.05)
        
    # Limits
    plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt_u).to_value('MeV cm-2 s-1'), color='blue',
             linewidth=2, linestyle='--', label=str(conf)+'% C.L.')
    plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt_d).to_value('MeV cm-2 s-1'), color='blue',
             linewidth=2, linestyle='--')
    plt.fill_between(E.to_value('GeV'), (E**2*dN_dEdSdt_u).to_value('MeV cm-2 s-1'),
                     (E**2*dN_dEdSdt_d).to_value('MeV cm-2 s-1'), color='blue', alpha=0.2)
    
    plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt_u).to_value('MeV cm-2 s-1'),
             color='magenta', linewidth=2, linestyle='--')
    plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt_d).to_value('MeV cm-2 s-1'),
             color='magenta', linewidth=2, linestyle='--')
    plt.fill_between(E.to_value('GeV'), (E**2*dNIC_dEdSdt_u).to_value('MeV cm-2 s-1'),
                     (E**2*dNIC_dEdSdt_d).to_value('MeV cm-2 s-1'), color='magenta', alpha=0.2)
    
    # Best fit
    plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt).to_value('MeV cm-2 s-1'),
             color='blue', linewidth=3, label='Best-fit model (Hadronic)')
    plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt).to_value('MeV cm-2 s-1'),
             color='magenta', linewidth=3, linestyle='-',label='Best-fit model (IC)')

    # Base
    plt.plot(E.to_value('GeV'), (E**2*base_dN_dEdSdt).to_value('MeV cm-2 s-1'),
             color='green', linewidth=3, label='Base model (Hadronic)')
    plt.plot(E.to_value('GeV'), (E**2*base_dNIC_dEdSdt).to_value('MeV cm-2 s-1'),
             color='chartreuse', linewidth=3, linestyle='-',label='Base model (IC)')

    plt.fill_between([30, 100e3], [0,0], [1e6,1e6], color='green', alpha=0.1, label='CTA energy range')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (GeV)')
    plt.ylabel(r'$\frac{E^2 dN}{dEdSdt}$ (MeV cm$^{-2}$ s$^{-1}$)')
    plt.xlim(1e-2, 2e5)
    plt.ylim(1e-10, 5e-6)
    plt.legend(fontsize=12)
    plt.savefig(cluster.output_dir+'/'+model_case+'_Gamma_spectrum.pdf')
    plt.close()

    #----- Profile
    fig = plt.figure(0, figsize=(8, 6))
    # MC
    for imc in range(Nmc):
        if imc == 0:
            plt.plot(r.to_value('kpc'), (prof_g_mc[imc]).to_value('cm-2 s-1 deg-2'),
                     color='blue', alpha=0.05, label='Monte Carlo')
        else:
            plt.plot(r.to_value('kpc'), (prof_g_mc[imc]).to_value('cm-2 s-1 deg-2'),
                     color='blue', alpha=0.05)
    for imc in range(Nmc):
        plt.plot(r.to_value('kpc'), (prof_ic_mc[imc]).to_value('cm-2 s-1 deg-2'),
                 color='magenta', alpha=0.05)

    # Limits
    plt.plot(r.to_value('kpc'), (dN_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
             color='blue', linewidth=2, linestyle='--', label=str(conf)+'% C.L.')
    plt.plot(r.to_value('kpc'), (dN_dSdtdO_d).to_value('cm-2 s-1 deg-2'),
             color='blue', linewidth=2, linestyle='--')
    plt.fill_between(r.to_value('kpc'), (dN_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
                     (dN_dSdtdO_d).to_value('cm-2 s-1 deg-2'), color='blue', alpha=0.2)
    
    plt.plot(r.to_value('kpc'), (dNIC_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
             color='magenta', linewidth=2, linestyle='--')
    plt.plot(r.to_value('kpc'), (dNIC_dSdtdO_d).to_value('cm-2 s-1 deg-2'),
             color='magenta', linewidth=2, linestyle='--')
    plt.fill_between(r.to_value('kpc'), (dNIC_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
                     (dNIC_dSdtdO_d).to_value('cm-2 s-1 deg-2'), color='magenta', alpha=0.2)

    # Best-fit
    plt.plot(r.to_value('kpc'), (dN_dSdtdO).to_value('cm-2 s-1 deg-2'),
             color='blue', linewidth=3, label='Best-fit model (Hadronic)')
    plt.plot(r.to_value('kpc'), (dNIC_dSdtdO).to_value('cm-2 s-1 deg-2'),
             color='magenta', linewidth=3, linestyle='-', label='Best-fit model (IC)')

    # Base
    plt.plot(r.to_value('kpc'), (base_dN_dSdtdO).to_value('cm-2 s-1 deg-2'),
             color='green', linewidth=3, label='Base model (Hadronic)')
    plt.plot(r.to_value('kpc'), (base_dNIC_dSdtdO).to_value('cm-2 s-1 deg-2'),
             color='chartreuse', linewidth=3, linestyle='-', label='Base model (IC)')
    
    plt.vlines((0.05*u.deg*cluster.cosmo.kpc_proper_per_arcmin(cluster.redshift)).to_value('kpc'),
               0,1, linestyle=':', color='k', label='CTA PSF (1 TeV)')
    plt.vlines(cluster.R500.to_value('kpc'), 0,1, linestyle='--', color='k', label='$R_{500}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'$\frac{dN}{dSdtd\Omega}$ (cm$^{-2}$ s$^{-1}$ deg$^{-2}$)')
    plt.xlim(10,5e3)
    plt.ylim(1e-16,1e-9)
    plt.legend(fontsize=12)
    plt.savefig(cluster.output_dir+'/'+model_case+'_Gamma_profile.pdf')
    plt.close()

    #----- Flux
    file = open(cluster.output_dir+'/'+model_case+'_Gamma_flux.txt','w')
    file.write('Energy E_gamma > 150 GeV \n')
    file.write('Hadronic        - Best  '+str(fluxH)+' cm-2 s-1 \n')
    file.write('                  Upper '+str(fluxH_u)+' cm-2 s-1 \n')
    file.write('                  Lower '+str(fluxH_d)+' cm-2 s-1 \n')
    file.write('Inverse Compton - Best  '+str(fluxIC)+' cm-2 s-1 \n')
    file.write('                  Upper '+str(fluxIC_u)+' cm-2 s-1 \n')
    file.write('                  Lower '+str(fluxIC_d)+' cm-2 s-1 \n')
    file.close()

    plotting.seaborn_1d(np.array(flux_g_mc)*1e12,
                        output_fig=cluster.output_dir+'/'+model_case+'_Gamma_flux_H.pdf',
                        ci=0.68, truth=None, best=fluxH*1e12,
                        label='$F_{\gamma,\ hadronic}$ ($10^{-12}$ s$^{-1}$ cm$^{-2}$)',
                        gridsize=100, alpha=(0.2, 0.4), 
                        figsize=(10,10), fontsize=12,
                        cols=[('blue','grey', 'orange')])
    plt.close("all")

    plotting.seaborn_1d(np.array(flux_ic_mc)*1e12,
                        output_fig=cluster.output_dir+'/'+model_case+'_Gamma_flux_IC.pdf',
                        ci=0.68, truth=None, best=fluxIC*1e12,
                        label='$F_{\gamma,\ IC}$ ($10^{-12}$ s$^{-1}$ cm$^{-2}$)',
                        gridsize=100, alpha=(0.2, 0.4), 
                        figsize=(10,10), fontsize=12,
                        cols=[('blue','grey', 'orange')])
    plt.close("all")
    print('  --> Gamma plots done')


#========================================
# Leptonic model
#========================================

def model_leptonic(params, cluster, data):
    '''
    Compute the leptonic model

    Parameters
    ----------
    - params (list of float): model parameters to be tested
    - cluster (minot object): cluster object
    - data (dict): radio data

    Output
    ------
    - p_synch (quantity array): radio profile
    - s_synch (quantity array): radio spectrum
    - i_synch (quantity array): radio spectral index profile

    '''

    #--- Extract parameters
    CR_X_E   = params[0]*1e-5
    CR_eta   = params[1]
    CR_index = params[2]
    Norm     = params[3]

    #--- Set parameters
    cluster = perseus_model_library.set_pure_leptonic_model(cluster, ('density', CR_eta), CR_X_E, CR_index)
    if app_steady: cluster.cre1_loss_model = 'Steady'

    #--- Profile
    r_synch, p_synch = cluster.get_synchrotron_profile(data['profile']['radius'], freq0=data['info']['prof_freq'])
    
    #--- Spectrum
    sfreq_bis = np.append(np.array([1]), data['spectrum']['freq'].to_value('MHz'))*u.MHz
    s_synch = cluster.get_synchrotron_spectrum(sfreq_bis, Rmin=data['info']['spec_Rmin'],
                                               Rmax=data['info']['spec_Rmax'], type_integral='cylindrical')[1]
    s_synch = s_synch[1:]
    
    #--- Spectral index
    if fit_index:
        p_synch1 = cluster.get_synchrotron_profile(data['index']['radius'], freq0=data['info']['idx_freq1'])[1]
        p_synch2 = cluster.get_synchrotron_profile(data['index']['radius'], freq0=data['info']['idx_freq2'])[1]
        upper = np.log10((p_synch1/p_synch2).to_value(''))
        lower = np.log10((data['info']['idx_freq1']/data['info']['idx_freq2']).to_value(''))
        i_synch  = -upper / lower
    else:
        i_synch = np.nan
        
    return Norm*p_synch, s_synch, i_synch

    
#========================================
# Hadronic model
#========================================

def model_hadronic(params, cluster, data):
    '''
    Compute the hadronic model

    Parameters
    ----------
    - params (list of float): model parameters to be tested
    - cluster (minot object): cluster object
    - data (dict): radio data

    Output
    ------
    - p_synch (quantity array): radio profile
    - s_synch (quantity array): radio spectrum
    - i_synch (quantity array): radio spectral index profile

    '''

    #--- Extract parameters
    CR_X_E   = params[0]*1e-2
    CR_eta   = params[1]
    CR_slope = params[2]
    Norm     = params[3]
    
    #--- Set parameters
    cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', CR_eta), CR_X_E, CR_slope)
    
    #--- Profile
    r_synch, p_synch = cluster.get_synchrotron_profile(data['profile']['radius'], freq0=data['info']['prof_freq'])
    
    #--- Spectrum
    sfreq_bis = np.append(np.array([1]), data['spectrum']['freq'].to_value('MHz'))*u.MHz
    s_synch = cluster.get_synchrotron_spectrum(sfreq_bis, 
                                               Rmin=data['info']['spec_Rmin'], Rmax=data['info']['spec_Rmax'], 
                                               type_integral='cylindrical')[1]
    s_synch = s_synch[1:]
    
    #--- Spectral index
    if fit_index:
        p_synch1 = cluster.get_synchrotron_profile(data['index']['radius'], freq0=data['info']['idx_freq1'])[1]
        p_synch2 = cluster.get_synchrotron_profile(data['index']['radius'], freq0=data['info']['idx_freq2'])[1]
        upper = np.log10((p_synch1/p_synch2).to_value(''))
        lower = np.log10((data['info']['idx_freq1']/data['info']['idx_freq2']).to_value(''))
        i_synch  = -upper/lower
    else:
        i_synch = np.nan
    
    return Norm*p_synch, s_synch, i_synch


#==================================================
# MCMC: Defines log prior
#==================================================

def lnprior(params, par_min, par_max):
    '''
    Return the flat prior on parameters

    Parameters
    ----------
    - params (list): the parameters
    - par_min (list): the minimum value for params
    - par_max (list): the maximum value for params

    Output
    ------
    - prior (float): the value of the prior, either 0 or -inf

    '''

    prior = 0.0
    
    for i in range(len(params)):
        if params[i] <= par_min[i] or params[i] >= par_max[i] :
            prior = -np.inf
            
    return prior


#==================================================
# MCMC: Defines log gaussian prior
#==================================================

def lngprior(params, par_gprior):
    '''
    Return the gaussian prior on parameters

    Parameters
    ----------
    - params (list): the parameters
    - par_gprior (list): the mean and sigma value for params

    Output
    ------
    - prior (float): the value of the prior

    '''
    
    prior = 0.0
    
    for i in range(len(params)):
        if par_gprior[1][i] != np.inf:
            norm = np.log(1.0/(np.sqrt(2*np.pi)*par_gprior[1][i]))
            prior += norm - 0.5*(params[i]-par_gprior[0][i])**2 / par_gprior[1][i]**2

    return prior


#========================================
# Hadronic log likelihood function
#========================================

def lnlike(params, cluster, data, par_min, par_max, par_gprior,
           fit_index=False, model_case='Hadronic'):
    '''
    log likelihood function

    Parameters
    ----------
    - params (array): grid of models
    - cluster (list): 

    Output
    ------
    - log likelihood value
    '''
    
    #--- Flat priors
    prior = lnprior(params, par_min, par_max)
    if prior == -np.inf: # Should not go for the model if prior out
        return -np.inf
    if np.isinf(prior):
        return -np.inf
    if np.isnan(prior):
        return -np.inf

    #--- Gaussian priors
    gprior = lngprior(params, par_gprior)
    
    #--- Model
    if model_case == 'Hadronic':
        prof_mod, spec_mod, idx_mod = model_hadronic(params, cluster, data)
    if model_case == 'Leptonic':
        prof_mod, spec_mod, idx_mod = model_leptonic(params, cluster, data)

    #--- Chi2
    # Profile chi2
    wgood1 = (data['profile']['radius']>data['info']['prof_Rmin'])
    wgood2 = (data['profile']['radius']<data['info']['prof_Rmax'])
    wgood = wgood1 * wgood2
    prof_resid = (data['profile']['flux']-prof_mod).to_value('Jy arcmin-2')[wgood]
    prof_chi2 = prof_resid**2/data['profile']['error'].to_value('Jy arcmin-2')[wgood]**2
    
    # Spectrum chi2
    spec_resid = (data['spectrum']['flux'].to_value('Jy')-spec_mod.to_value('Jy'))**2
    spec_chi2 = spec_resid/data['spectrum']['error'].to_value('Jy')**2
    
    # Spectral index chi2
    if fit_index:
        wgood1 = (data['index']['radius']>data['info']['idx_Rmin'])
        wgood2 = (data['index']['radius']<data['info']['idx_Rmax'])
        wgood = wgood1 * wgood2
        idx_chi2 = ((data['index']['idx'][wgood]-idx_mod[wgood])**2)/data['index']['error'][wgood]**2
    else:
        idx_chi2 = np.nan
        
    # Chi2 tot
    chi2_tot = -0.5*np.nansum(prof_chi2) - 0.5*np.nansum(spec_chi2)
    if fit_index: 
        chi2_tot = chi2_tot - 0.5*np.nansum(idx_chi2)

    return chi2_tot + gprior


#========================================
# Standard curve fit first
#========================================

def run_curvefit(cluster, radio_data, par0, par_min, par_max,
                 fit_index=False, model_case='Hadronic'):
    '''
    Run the curvefit

    Parameters
    ----------
    - cluster (minot object): cluster object
    - par0 (list): guess parameters
    - par_min/max (list): flat prior limit on parameters
    - mcmc_nsteps (int): number of MCMC steps
    - nwalkers (int): number of walkers
    - run_mcmc (bool): run the MCMC or not
    - reset_mcmc (bool): reset MCMC or start from existing chains
    - fit_index (bool): fit the spectral index
    - model_case (string): case considered, 'Hadronic' or 'Leptonic'

    Output
    ------
    - The optimal curvefit parameters
    '''

    # fit function
    if len(par0 == 3):
        def fitfunc(xdata, p1,p2,p3):
            params = [p1,p2,p3]
            if model_case == 'Hadronic':
                prof_mod, spec_mod, idx_mod = model_hadronic(params, cluster, radio_data)
            if model_case == 'Leptonic':
                prof_mod, spec_mod, idx_mod = model_leptonic(params, cluster, radio_data)
            model = np.append(prof_mod.to_value('Jy arcmin-2'), spec_mod.to_value('Jy'))
            if fit_index:
                model = np.append(model,idx_mod)
            return model
        
    if len(par0) == 4:
        def fitfunc(xdata, p1,p2,p3,p4):
            params = [p1,p2,p3,p4]
            if model_case == 'Hadronic':
                prof_mod, spec_mod, idx_mod = model_hadronic(params, cluster, radio_data)
            if model_case == 'Leptonic':
                prof_mod, spec_mod, idx_mod = model_leptonic(params, cluster, radio_data)
            model = np.append(prof_mod.to_value('Jy arcmin-2'), spec_mod.to_value('Jy'))
            if fit_index:
                model = np.append(model,idx_mod)
            return model

    if len(par0) == 5:
        def fitfunc(xdata, p1,p2,p3,p4):
            params = [p1,p2,p3,p4]
            if model_case == 'Hadronic':
                prof_mod, spec_mod, idx_mod = model_hadronic(params, cluster, radio_data)
            if model_case == 'Leptonic':
                prof_mod, spec_mod, idx_mod = model_leptonic(params, cluster, radio_data)
            model = np.append(prof_mod.to_value('Jy arcmin-2'), spec_mod.to_value('Jy'))
            if fit_index:
                model = np.append(model,idx_mod)
            return model
        
    # Define the data
    Np = len(radio_data['profile']['radius'].to_value('kpc'))
    Ns = len(radio_data['spectrum']['freq'].to_value('MHz'))
    Ni = len(radio_data['index']['idx'])
    xdata = np.append(radio_data['profile']['radius'].to_value('kpc'),
                      radio_data['spectrum']['freq'].to_value('MHz'))
    ydata = np.append(radio_data['profile']['flux'].to_value('Jy arcmin-2'),
                      radio_data['spectrum']['flux'].to_value('Jy'))
    sigma = np.append(radio_data['profile']['error'].to_value('Jy arcmin-2'),
                      radio_data['spectrum']['error'].to_value('Jy'))
    #sigma = np.zeros(Np+Ns)+1.0
    wbad1 = (radio_data['profile']['radius']<radio_data['info']['prof_Rmin'])
    wbad2 = (radio_data['profile']['radius']>radio_data['info']['prof_Rmax'])
    wbad3 = (radio_data['spectrum']['freq'].to_value('MHz') < 0)
    wbad = np.append(wbad1*wbad2, wbad3)
    if fit_index:
        xdata = np.append(xdata, radio_data['index']['radius'].to_value('kpc'))
        ydata = np.append(ydata, radio_data['index']['idx'])
        sigma = np.append(sigma, radio_data['index']['error'])
        #sigma = np.zeros(Np+Ns+Ni)+1.0
        wbad4 = (radio_data['index']['radius']<radio_data['info']['idx_Rmin'])
        wbad5 = (radio_data['index']['radius']>radio_data['info']['idx_Rmax'])
        wbad = np.append(wbad, wbad4*wbad5)
    sigma[wbad] = sigma[wbad] * 1e10
    
    # fit
    p_opt, p_cov = curve_fit(fitfunc, xdata, ydata,
                             p0=par0, sigma=sigma, absolute_sigma=True,
                             check_finite=True, bounds=(par_min, par_max),
                             method=None, jac=None)
    
    print('       parameters:')
    print(p_opt)
    print('       covariance^1/2:')
    print(p_cov**0.5)
    print('')

    return p_opt, p_cov


#========================================
# Run the MCMC
#========================================

def run_function_mcmc(cluster, radio_data, par0, par_min, par_max, par_gprior,
                      mcmc_nsteps=1000, nwalkers=10, moves=None,
                      run_mcmc=True, reset_mcmc=False,
                      fit_index=False, model_case='Hadronic'):
    '''
    Run the MCMC

    Parameters
    ----------
    - cluster (minot object): cluster object
    - par0 (list): guess parameters
    - par_min/max (list): flat prior limit on parameters
    - mcmc_nsteps (int): number of MCMC steps
    - nwalkers (int): number of walkers
    - run_mcmc (bool): run the MCMC or not
    - reset_mcmc (bool): reset MCMC or start from existing chains
    - fit_index (bool): fit the spectral index
    - model_case (string): case considered, 'Hadronic' or 'Leptonic'

    Output
    ------
    - The chains are produced and saved
    '''

    
    #---------- Check if a MCMC sampler was already recorded
    sampler_exist = os.path.exists(cluster.output_dir+'/'+model_case+'_sampler.h5')
    if sampler_exist:
        print('    Existing sampler: '+cluster.output_dir+'/'+model_case+'_sampler.h5')
    else:
        print('    No existing sampler found')
    
    #---------- MCMC parameters
    ndim = len(par0)
    if nwalkers < 2*ndim:
        print('    nwalkers should be at least twice the number of parameters.')
        print('    nwalkers --> '+str(ndim*2))
        print('')
        nwalkers = ndim*2
            
    #----- Define the MCMC
    backend = emcee.backends.HDFBackend(cluster.output_dir+'/'+model_case+'_sampler.h5')
    pos = mcmc_common.chains_starting_point(par0, 0.1, par_min, par_max, nwalkers)
    if sampler_exist:
        if reset_mcmc:
            print('    Reset MCMC even though sampler already exists')
            backend.reset(nwalkers, ndim)
        else:
            print('    Use existing MCMC sampler')
            print("    --> Initial size: {0}".format(backend.iteration))
            pos = None
    else:
        print('    No pre-existing sampler, start from scratch')
        backend.reset(nwalkers, ndim)
        
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike,
                                    args=[cluster, radio_data, par_min, par_max, par_gprior,
                                          fit_index, model_case],
                                    pool=Pool(cpu_count()), moves=moves,
                                    backend=backend)

    #----- Run the MCMC
    if run_mcmc:
        print('    Runing '+str(mcmc_nsteps)+' MCMC steps')
        res = sampler.run_mcmc(pos, mcmc_nsteps, progress=True)
    
        #----- Save sampler
        with open(cluster.output_dir+'/'+model_case+'_sampler.pkl', 'wb') as output:
            pickle.dump(sampler, output, pickle.HIGHEST_PROTOCOL)
        
        
#========================================
# Main function
#========================================

if __name__ == "__main__":
    
    #========== Parameters
    Nmc         = 100              # Number of Monte Carlo trials
    fit_index   = False            # Fit the spectral index profile
    app_steady  = True             # Application of steady state losses
    mcmc_nsteps = 89             # number of MCMC points
    mcmc_nwalk  = 96#2*cpu_count()    # number of walkers
    mcmc_burnin = 2000             # number of MCMC burnin points
    mcmc_reset  = False            # Reset the MCMC
    run_mcmc    = False             # Run the MCMC
    basedata    = 'Pedlar1990'     # 'Gitti2002', 'Pedlar1990'
    model_case  = 'Hadronic'       # 'Hadronic' or 'Leptonic'
    mag_case    = 'Taylor2006'
    #mag_case    = 'Bonafede2010up'
    #mag_case    = 'Walker2017'
    output_dir = '/sps/cta/llr/radam/PerseusGammaCalib'
    #output_dir  = '/Users/adam/Project/CTA/Phys/Outputs/Perseus_KSP_calibration/Calib'
    output_dir = output_dir+'_'+model_case+'_'+mag_case+'_'+basedata

    #========== Information
    print('========================================')
    print(' Starting the CR calibration of Perseus ')
    print('========================================')
    print('   Number of CPU available:      ', cpu_count())
    print('   Number of MC planed:          ', Nmc)
    print('   Burnin:                       ', mcmc_burnin)
    print('   Reset MCMC?:                  ', mcmc_reset)
    print('   Run MCMC?:                    ', run_mcmc)
    print('   Number of walkers:            ', mcmc_nwalk)
    print('   Number of steps:              ', mcmc_nsteps)
    print('   Fitting for spectral index?:  ', fit_index)
    print('   Applying SS loss (leptonic)?: ', app_steady)
    print('----------------------------------------')
    print('')

    #========== Make directory
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    print('-----> Working in '+output_dir)
    
    #========== Define the cluster model
    if model_case == 'Hadronic':
        cluster = perseus_model_library.default_model(directory=output_dir)
        cluster = perseus_model_library.set_magnetic_field_model(cluster, case=mag_case)
        cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', 1.0), 1e-2, 2.5)
        cluster.Npt_per_decade_integ = 10
    elif model_case == 'Leptonic':
        cluster = perseus_model_library.default_model(directory=output_dir)
        cluster = perseus_model_library.set_magnetic_field_model(cluster, case=mag_case)
        cluster = perseus_model_library.set_pure_leptonic_model(cluster, ('density', 1.0), 1e-5, 2.0)
        if app_steady: cluster.cre1_loss_model = 'Steady'
        cluster.Npt_per_decade_integ = 30
    else:
        raise ValueError('Only Hadronic or Leptonic are possible')
    
    #========== Data
    print('')
    print('-----> Getting the radio data')
    radio_data = perseus_data_library.get_radio_data(cluster.cosmo, cluster.redshift, prof_file=basedata)
    
    #========== Starting point
    print('')
    print('-----> Get guess parameter with curvefit')
    if model_case == 'Hadronic':
        param_name = ['X_{CRp} (x10^{-2})', '\\eta_{CRp}', '\\alpha_{CRp}', 'Norm']
        par0       = [1.0,     1.0,    2.5,    1.0]
        par_min    = [0.0,     0.0,    2.0,    0.5]
        par_max    = [20,      5.0,    4.0,    1.5]
        par_gprior = ([0.0,    1.0,    2.5,    1.0],
                      [20., np.inf, np.inf,    0.1])
    if model_case == 'Leptonic':
        param_name = ['X_{CRe} (x10^{-5})', '\\eta_{CRe}', '\\alpha_{CRe}', 'Norm']
        par0       = [1.0,     1.0,    2.0,    1.0]
        par_min    = [0,       0,      2,      0.5]
        par_max    = [1e4,     5,      5,      1.5]
        par_gprior = ([1.0,    1.0,    2.5,    1.0],
                      [np.inf, np.inf, np.inf, 0.1])
        
    if mcmc_reset and run_mcmc:
        par_opt, par_cov = run_curvefit(cluster, radio_data, np.array(par0), par_min, par_max,
                                        fit_index=fit_index, model_case=model_case)
    else:
        par_opt = np.array(par0) # this is not used in this case, but just to give something
        # Get the covariance
        with open(cluster.output_dir+'/'+model_case+'_sampler.pkl', 'rb') as f:
            sampler = pickle.load(f)
        param_chains = sampler.chain[:, mcmc_burnin:, :]
        chain_list = []
        for i in range(param_chains.shape[2]):
            chain_list.append(param_chains[:,:,i].flatten())
        par_cov = np.cov(np.array(chain_list))
        print('---- Param chain covariance^1/2:')
        print(par_cov**0.5)
        print('')

    #========== MCMC fit
    print('')
    print('-----> Going for the MCMC fit')
    moves = emcee.moves.StretchMove(a=2.0)
    run_function_mcmc(cluster, radio_data, par_opt, par_min, par_max, par_gprior,
                      mcmc_nsteps=mcmc_nsteps, nwalkers=mcmc_nwalk, moves=moves,
                      run_mcmc=run_mcmc, reset_mcmc=mcmc_reset,
                      fit_index=fit_index, model_case=model_case)
    
    #========== Post analysis
    print('')
    print('-----> Starting the post-analysis')
    post_analysis(cluster, radio_data, param_name,
                  par_min, par_max, mcmc_burnin,
                  conf=68.0, Nmc=Nmc, model_case=model_case)
    
