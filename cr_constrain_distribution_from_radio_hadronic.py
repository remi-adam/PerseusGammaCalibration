"""
This script is part of the CTA key science project dedicated to Perseus. 
Here, we constrain the Perseus CR distribution (spectral and spatial) using 
radio data, via an MCMC approach. This is the pure hadronic model.
"""

#========================================
# Imports
#========================================

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

from PerseusGammaCalibration import perseus_model_library
from PerseusGammaCalibration import perseus_data_library
from kesacco.Tools.plotting import seaborn_corner
import minot

# Modify plotting parameters
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
# Functions
#========================================

#========== Defines log likelihood
def lnlike(params, cluster, data):
    #--- Priors
    cond1 = params[0]>=0  and params[0]<=100  # X
    cond2 = params[1]>=-5 and params[1]<=5    # eta_CRp
    cond3 = params[2]>=2  and params[2]<=4    # Slope
    if cond1 and cond2 and cond3:
        
        prof_mod, spec_mod, idx_mod = model(params, cluster, data)
        
        # Profile chi2
        wgood = (data['profile']['radius']>data['info']['prof_Rmin'])*(data['profile']['radius']<data['info']['prof_Rmax'])
        prof_resid = (data['profile']['flux'].to_value('Jy arcmin-2')[wgood]-prof_mod.to_value('Jy arcmin-2')[wgood])**2
        prof_chi2 = prof_resid/data['profile']['error'].to_value('Jy arcmin-2')[wgood]**2

        # Spectrum chi2
        spec_resid = (data['spectrum']['flux'].to_value('Jy')-spec_mod.to_value('Jy'))**2
        spec_chi2 = spec_resid/data['spectrum']['error'].to_value('Jy')**2

        # Spectral index chi2
        wgood = (data['index']['radius']>data['info']['idx_Rmin'])*(data['index']['radius']<data['info']['idx_Rmax'])
        idx_chi2 = ((data['index']['idx'][wgood]-idx_mod[wgood])**2)/data['index']['error'][wgood]**2
        
        # Chi2 tot
        chi2_tot = -0.5*np.nansum(prof_chi2) - 0.5*np.nansum(spec_chi2)
        if fit_index: 
            chi2_tot = chi2_tot - 0.5*np.nansum(idx_chi2)
        
        return chi2_tot
    else:
        #print('Prior conditions are ', cond1, cond2, cond3, cond4)
        return -np.inf

#========== Defines model
def model(params, cluster, data):
    #--- Extract parameters
    CR_X_E   = params[0]*1e-2
    CR_eta   = params[1]
    CR_slope = params[2]
    
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
    p_synch1 = cluster.get_synchrotron_profile(data['index']['radius'], freq0=data['info']['idx_freq1'])[1]
    p_synch2 = cluster.get_synchrotron_profile(data['index']['radius'], freq0=data['info']['idx_freq2'])[1]
    upper = np.log10((p_synch1/p_synch2).to_value(''))
    lower = np.log10((data['info']['idx_freq1']/data['info']['idx_freq2']).to_value(''))
    i_synch  = -upper/lower
    
    return p_synch, s_synch, i_synch


#========================================
# Main function
#========================================

if __name__ == "__main__":

    #========== Parameters
    Nmc         = 10         # Number of Monte Carlo trials
    fit_index   = False       # Fit the spectral index profile
    mcmc_nsteps = 100        # number of MCMC points
    mcmc_burnin = 10         # number of MCMC burnin points
    mcmc_reset  = False      # Reset the MCMC
    
    #========== Define the cluster model
    cluster = perseus_model_library.default_model()
    cluster = perseus_model_library.set_magnetic_field_model(cluster, case='Taylor2006')
    cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', 2.0), 3e-3, 2.5)
    cluster.Npt_per_decade_integ = 10
    
    #========== Data
    radio_data = perseus_data_library.get_radio_data(cluster.cosmo, cluster.redshift)

    #========== MCMC fit
    #----- Define the MCMC
    par0 = np.array([0.3, 2.0, 2.5])
    
    param_name = [r'$X_{CRp}$(%)', r'$\eta_{CRp}$', r'$\alpha_{CRp}$']
    ndim, nwalkers, nsteps, burnin = len(par0), 10, mcmc_nsteps, mcmc_burnin
    pos = [par0 + par0*1e-1*np.random.randn(ndim) for i in range(nwalkers)]
    
    #----- Define the MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, args=[cluster, radio_data], pool=Pool(cpu_count()))
    
    #----- Restart from where it is or reset
    if reset:
        sampler.reset()
    
    #----- Rune the MCMC
    res = sampler.run_mcmc(pos, nsteps, progress=True)
    
    #----- Save sampler
    with open(cluster.output_dir+'/Hadronic_sampler.pkl', 'wb') as output:
        pickle.dump(sampler, output, pickle.HIGHEST_PROTOCOL)

    #----- Burn in
    param_chains = sampler.chain[:, burnin:, :]
    lnL_chains = sampler.lnprobability[:, burnin:]
    
    #----- Get the best fit parameters
    wbest = (lnL_chains == np.amax(lnL_chains))
    param_best = []
    for i in range(ndim):
        param_best.append(((param_chains[:,:,i])[wbest])[0])
        
    #----- MC parameters
    param_flat = param_chains.reshape(param_chains.shape[0]*param_chains.shape[1],param_chains.shape[2])
    Nsample = len(param_flat[:,0])-1
    param_MC = np.zeros((Nmc, ndim))
    for i in range(Nmc):
        param_MC[i,:] = param_flat[np.random.randint(0, high=Nsample), :] # randomly taken from chains

    #========== Results
    #----- Chains
    for i in range(ndim):
        plt.figure(i)
        plt.plot((param_chains[:,:,i]).flatten())
        plt.xlabel('Step number')
        plt.ylabel(param_name[i])
        
    #----- Parameters
    print(param_best)

    #----- Parameter space
    par_flat = param_chains.reshape(param_chains.shape[0]*param_chains.shape[1], param_chains.shape[2])
    fig = corner.corner(par_flat, bins=50, color='k', smooth=2, 
                        labels=param_name, quantiles=(0.16, 0.84), use_math_text=True)
    fig.savefig(cluster.output_dir+'/Hadronic_triangle_corner.pdf')
    
    df = pd.DataFrame(par_flat, columns=param_name)
    seaborn_corner(df, n_levels=30, cols=[('royalblue', 'k', 'grey', 'Blues')], 
                   ci2d=[0.68, 0.95], gridsize=100, zoom=0.1,
                   linewidth=2.0, alpha=(0.3, 0.3, 1.0), figsize=(9,9))
    
    #----- Data versus model
    # Best-fit
    prof_best, spec_best, idx_best = model(param_best, cluster, radio_data)
    
    # MC sampling
    prof_mc = []
    spec_mc = []
    idx_mc  = []
    
    for imc in range(Nmc):
        prof_mci, spec_mci, idx_mci = model(param_MC[imc,:], cluster, radio_data)
        prof_mc.append(prof_mci)
        spec_mc.append(spec_mci)
        idx_mc.append(idx_mci)
        
    # Limits
    prof_u = np.percentile(np.array(prof_mc), 100-(100-68)/2.0, axis=0)*prof_mc[0].unit
    prof_d = np.percentile(np.array(prof_mc), (100-68)/2.0, axis=0)*prof_mc[0].unit
    
    spec_u = np.percentile(np.array(spec_mc), 100-(100-68)/2.0, axis=0)*spec_mc[0].unit
    spec_d = np.percentile(np.array(spec_mc), (100-68)/2.0, axis=0)*spec_mc[0].unit
    
    idx_u = np.percentile(np.array(idx_mc), 100-(100-68)/2.0, axis=0)
    idx_d = np.percentile(np.array(idx_mc), (100-68)/2.0, axis=0)
    
    
    fig = plt.figure(0, figsize=(18, 16))

    #----- Spectrum
    ax = plt.subplot(221)
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_mc[imc].to_value('Jy'),
                 color='blue', alpha=0.05)
        
    # Limits
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_u.to_value('Jy'), color='blue',
             linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_d.to_value('Jy'), color='blue',
             linewidth=2, linestyle='--')
    plt.fill_between(radio_data['spectrum']['freq'].to_value('MHz'), spec_d.to_value('Jy'), spec_u.to_value('Jy'),
                     color='blue', alpha=0.2)
    
    # Best fit and data
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_best.to_value('Jy'), color='blue',
             linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['spectrum']['freq'].to_value('MHz'), radio_data['spectrum']['flux'].to_value('Jy'),
                 radio_data['spectrum']['error'].to_value('Jy'),
                 marker='o', color='red', linestyle='', label='Data')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('frequency (MHz)')
    plt.ylabel('flux (Jy)')
    plt.legend()
    plt.grid()
    
    #----- Profile
    ax = plt.subplot(222)
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_mc[imc].to_value('Jy arcmin-2'),
                 color='blue', alpha=0.05)

    # Limits    
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_u.to_value('Jy arcmin-2'),
             color='blue', linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_d.to_value('Jy arcmin-2'),
             color='blue', linewidth=2, linestyle='--')
    plt.fill_between(radio_data['profile']['radius'].to_value('kpc'), prof_d.to_value('Jy arcmin-2'),
                     prof_u.to_value('Jy arcmin-2'), color='blue', alpha=0.2)
    
    # Best and data
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_best.to_value('Jy arcmin-2'),
             color='blue', linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['profile']['radius'].to_value('kpc'), radio_data['profile']['flux'].to_value('Jy arcmin-2'), 
                 yerr=radio_data['profile']['error'].to_value('Jy arcmin-2'),
                 marker='o', linestyle='', color='red', label='Data')
    
    plt.plot([radio_data['info']['prof_Rmin'].to_value('kpc'),radio_data['info']['prof_Rmin'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.plot([radio_data['info']['prof_Rmax'].to_value('kpc'),radio_data['info']['prof_Rmax'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.fill_between([0,radio_data['info']['prof_Rmin'].to_value('kpc')], [0,0], [1e6,1e6], color='grey',
                     alpha=0.2, label='Excluded region')
    plt.fill_between([radio_data['info']['prof_Rmax'].to_value('kpc'), np.inf], [0,0], [1e6,1e6],
                     color='grey', alpha=0.2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('radius (kpc)')
    plt.ylabel('surface brightness (Jy/arcmin$^2$)')
    plt.xlim(10,250)
    plt.ylim(1e-3,1e1)
    plt.legend()
    plt.grid()
    
    #----- Spectral index
    ax = plt.subplot(224)
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_mc[imc],
                 color='blue', alpha=0.05)
        
    # Limits
    plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_u, color='blue',
             linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_d, color='blue',
             linewidth=2, linestyle='--')
    plt.fill_between(radio_data['index']['radius'].to_value('kpc'), idx_d, idx_u,
                     color='blue', alpha=0.2)
    
    # Best and data    
    plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_best, color='blue',
             linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['index']['radius'].to_value('kpc'), radio_data['index']['idx'],
                 yerr=radio_data['index']['error'],
                 marker='o', linestyle='', color='red', label='Data')
    
    plt.plot([radio_data['info']['idx_Rmin'].to_value('kpc'),radio_data['info']['idx_Rmin'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.plot([radio_data['info']['idx_Rmax'].to_value('kpc'),radio_data['info']['idx_Rmax'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.fill_between([0,radio_data['info']['idx_Rmin'].to_value('kpc')], [0,0], [1e6,1e6],
                     color='grey', alpha=0.1, label='Excluded region')
    plt.fill_between([radio_data['info']['idx_Rmax'].to_value('kpc'), np.inf],
                     [0,0], [1e6,1e6], color='grey', alpha=0.1)
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlabel('radius (kpc)')
    plt.ylabel('spectral index')
    plt.xlim(10,250)
    plt.ylim(0.5,2.5)
    plt.legend()
    plt.grid()
    
    #----- Spectrum
    fig = plt.figure(0, figsize=(8, 6))
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_mc[imc].to_value('Jy'),
                 color='blue', alpha=0.05)
        
    # Limits
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_u.to_value('Jy'),
             color='blue', linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_d.to_value('Jy'),
             color='blue', linewidth=2, linestyle='--')
    plt.fill_between(radio_data['spectrum']['freq'].to_value('MHz'), spec_d.to_value('Jy'),
                     spec_u.to_value('Jy'), color='blue', alpha=0.2)
    
    # Best fit and data
    plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_best.to_value('Jy'),
             color='blue', linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['spectrum']['freq'].to_value('MHz'), radio_data['spectrum']['flux'].to_value('Jy'),
                 radio_data['spectrum']['error'].to_value('Jy'),
                 marker='o', color='red', linestyle='', label='Data')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('frequency (MHz)')
    plt.ylabel('flux (Jy)')
    plt.legend()
    plt.savefig(cluster.output_dir+'/Haronic_Radio_spectrum.pdf')
    plt.close()
    
    #----- Profile
    fig = plt.figure(0, figsize=(8, 6))
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_mc[imc].to_value('Jy arcmin-2'),
                 color='blue', alpha=0.05)
        
    # Limits    
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_u.to_value('Jy arcmin-2'),
             color='blue', linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_d.to_value('Jy arcmin-2'),
             color='blue', linewidth=2, linestyle='--')
    plt.fill_between(radio_data['profile']['radius'].to_value('kpc'), prof_d.to_value('Jy arcmin-2'),
                     prof_u.to_value('Jy arcmin-2'), color='blue', alpha=0.2)
    
    
    # Best and data
    plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_best.to_value('Jy arcmin-2'),
             color='blue', linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['profile']['radius'].to_value('kpc'), radio_data['profile']['flux'].to_value('Jy arcmin-2'), 
                 yerr=radio_data['profile']['error'].to_value('Jy arcmin-2'),
                 marker='o', linestyle='', color='red', label='Data')
    
    plt.plot([radio_data['info']['prof_Rmin'].to_value('kpc'),radio_data['info']['prof_Rmin'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.plot([radio_data['info']['prof_Rmax'].to_value('kpc'),radio_data['info']['prof_Rmax'].to_value('kpc')],
             [0,1e6], linestyle='--', color='grey')
    plt.fill_between([0,radio_data['info']['prof_Rmin'].to_value('kpc')], [0,0], [1e6,1e6],
                     color='grey', alpha=0.2, label='Excluded region')
    plt.fill_between([radio_data['info']['prof_Rmax'].to_value('kpc'), np.inf], [0,0], [1e6,1e6],
                     color='grey', alpha=0.2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('radius (kpc)')
    plt.ylabel('surface brightness (Jy/arcmin$^2$)')
    plt.xlim(10,250)
    plt.ylim(1e-3,1e1)
    plt.savefig(cluster.output_dir+'/Haronic_Radio_profile.pdf')
    plt.close()
    
    #----- Spectral index
    fig = plt.figure(0, figsize=(8, 6))
    # MC
    for imc in range(Nmc):
        plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_mc[imc], color='blue', alpha=0.05)
        
    # Limits
    plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_u, color='blue',
             linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_d, color='blue',
             linewidth=2, linestyle='--')
    plt.fill_between(radio_data['index']['radius'].to_value('kpc'), idx_d, idx_u,
                     color='blue', alpha=0.2)
    
    # Best and data    
    plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_best, color='blue',
             linewidth=3, label='Best-fit')
    plt.errorbar(radio_data['index']['radius'].to_value('kpc'), radio_data['index']['idx'],
                 yerr=radio_data['index']['error'],
                 marker='o', linestyle='', color='red', label='Data')
    
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
    plt.savefig(cluster.output_dir+'/Haronic_Radio_index.pdf')
    plt.close()
    
    #========== Implication for gamma rays
    energy = np.logspace(-2,6,100)*u.GeV
    radius = np.logspace(0,4,100)*u.kpc
    
    # Best-fit
    cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', param_best[1]),
                                                            param_best[0]*1e-2, param_best[2])
    
    #----- Hadronic
    E, dN_dEdSdt = cluster.get_gamma_spectrum(energy, 
                                              Rmin=None, Rmax=cluster.R500,
                                              type_integral='cylindrical',
                                              Rmin_los=None, NR500_los=5.0)
    
    r, dN_dSdtdO = cluster.get_gamma_profile(radius, 
                                             Emin=50*u.GeV, Emax=100*u.TeV, 
                                             Energy_density=False, Rmin_los=None, NR500_los=5.0)
    
    #----- Inverse Compton
    E, dNIC_dEdSdt = cluster.get_ic_spectrum(energy, 
                                             Rmin=None, Rmax=cluster.R500,
                                             type_integral='cylindrical',
                                             Rmin_los=None, NR500_los=5.0)
    
    r, dNIC_dSdtdO = cluster.get_ic_profile(radius, 
                                            Emin=50*u.GeV, Emax=100*u.TeV, 
                                            Energy_density=False, Rmin_los=None, NR500_los=5.0)
    
    # Monte Carlo sampling: hadronic
    prof_g_mc = []
    spec_g_mc = []
    
    for imc in range(Nmc):
        cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', param_MC[imc,1]), 
                                                               param_MC[imc,0]*1e-2, param_MC[imc,2])
    
        spec_g_mci = cluster.get_gamma_spectrum(energy, Rmin=None, Rmax=cluster.R500,
                                                type_integral='cylindrical',
                                                Rmin_los=None, NR500_los=5.0)[1]
    
        prof_g_mci = cluster.get_gamma_profile(radius, Emin=50*u.GeV, Emax=100*u.TeV, 
                                               Energy_density=False, Rmin_los=None, NR500_los=5.0)[1]
        prof_g_mc.append(prof_g_mci)
        spec_g_mc.append(spec_g_mci)
        
        
        # Limits: hadronic
        dN_dEdSdt_u = np.percentile(np.array(spec_g_mc), 100-(100-68)/2.0, axis=0)*spec_g_mc[0].unit
        dN_dEdSdt_d = np.percentile(np.array(spec_g_mc), (100-68)/2.0, axis=0)*spec_g_mc[0].unit
        
        dN_dSdtdO_u = np.percentile(np.array(prof_g_mc), 100-(100-68)/2.0, axis=0)*prof_g_mc[0].unit
        dN_dSdtdO_d = np.percentile(np.array(prof_g_mc), (100-68)/2.0, axis=0)*prof_g_mc[0].unit
        
        # Monte Carlo sampling: IC
        prof_ic_mc = []
        spec_ic_mc = []
        
    for imc in range(Nmc):
        cluster = perseus_model_library.set_pure_hadronic_model(cluster, ('density', param_MC[imc,1]), 
                                                               param_MC[imc,0]*1e-2, param_MC[imc,2])
    
        spec_ic_mci = cluster.get_ic_spectrum(energy, Rmin=None, Rmax=cluster.R500,
                                              type_integral='cylindrical',
                                              Rmin_los=None, NR500_los=5.0)[1]
    
        prof_ic_mci = cluster.get_ic_profile(radius, Emin=50*u.GeV, Emax=100*u.TeV, 
                                             Energy_density=False, Rmin_los=None, NR500_los=5.0)[1]
        prof_ic_mc.append(prof_ic_mci)
        spec_ic_mc.append(spec_ic_mci)
    
    # Limits: IC
    dNIC_dEdSdt_u = np.percentile(np.array(spec_ic_mc), 100-(100-68)/2.0, axis=0)*spec_ic_mc[0].unit
    dNIC_dEdSdt_d = np.percentile(np.array(spec_ic_mc), (100-68)/2.0, axis=0)*spec_ic_mc[0].unit
    
    dNIC_dSdtdO_u = np.percentile(np.array(prof_ic_mc), 100-(100-68)/2.0, axis=0)*prof_ic_mc[0].unit
    dNIC_dSdtdO_d = np.percentile(np.array(prof_ic_mc), (100-68)/2.0, axis=0)*prof_ic_mc[0].unit
    
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
                 color='grey', alpha=0.05)
        
    # Limits
    plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt_u).to_value('MeV cm-2 s-1'), color='blue',
             linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt_d).to_value('MeV cm-2 s-1'), color='blue',
             linewidth=2, linestyle='--')
    plt.fill_between(E.to_value('GeV'), (E**2*dN_dEdSdt_u).to_value('MeV cm-2 s-1'),
                     (E**2*dN_dEdSdt_d).to_value('MeV cm-2 s-1'), color='blue', alpha=0.2)
    
    plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt_u).to_value('MeV cm-2 s-1'),
             color='k', linewidth=2, linestyle='--')
    plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt_d).to_value('MeV cm-2 s-1'),
             color='k', linewidth=2, linestyle='--')
    plt.fill_between(E.to_value('GeV'), (E**2*dNIC_dEdSdt_u).to_value('MeV cm-2 s-1'),
                     (E**2*dNIC_dEdSdt_d).to_value('MeV cm-2 s-1'), color='k', alpha=0.2)
    
    # Best fit
    plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt).to_value('MeV cm-2 s-1'),
             color='blue', linewidth=3, label='Best-fit model (Hadronic)')
    plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt).to_value('MeV cm-2 s-1'),
             color='k', linewidth=3, linestyle='-',label='Best-fit model (IC)')
    
    plt.fill_between([30, 100e3], [0,0], [1e6,1e6], color='red', alpha=0.1, label='CTA energy range')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (GeV)')
    plt.ylabel(r'$\frac{E^2 dN}{dEdSdt}$ (MeV cm$^{-2}$ s$^{-1}$)')
    plt.xlim(1e-2, 2e5)
    plt.ylim(1e-10, 5e-6)
    plt.legend(fontsize=14)
    plt.savefig(cluster.output_dir+'/Haronic_Gamma_spectrum.pdf')

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
                 color='grey', alpha=0.05)

    # Limits
    plt.plot(r.to_value('kpc'), (dN_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
             color='blue', linewidth=2, linestyle='--', label='68% C.L.')
    plt.plot(r.to_value('kpc'), (dN_dSdtdO_d).to_value('cm-2 s-1 deg-2'),
             color='blue', linewidth=2, linestyle='--')
    plt.fill_between(r.to_value('kpc'), (dN_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
                     (dN_dSdtdO_d).to_value('cm-2 s-1 deg-2'), color='blue', alpha=0.2)
    
    plt.plot(r.to_value('kpc'), (dNIC_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
             color='k', linewidth=2, linestyle='--')
    plt.plot(r.to_value('kpc'), (dNIC_dSdtdO_d).to_value('cm-2 s-1 deg-2'),
             color='k', linewidth=2, linestyle='--')
    plt.fill_between(r.to_value('kpc'), (dNIC_dSdtdO_u).to_value('cm-2 s-1 deg-2'),
                     (dNIC_dSdtdO_d).to_value('cm-2 s-1 deg-2'), color='k', alpha=0.2)

    # Best-fit
    plt.plot(r.to_value('kpc'), (dN_dSdtdO).to_value('cm-2 s-1 deg-2'),
             color='blue', linewidth=3, label='Best-fit model (Hadronic)')
    plt.plot(r.to_value('kpc'), (dNIC_dSdtdO).to_value('cm-2 s-1 deg-2'),
             color='k', linewidth=3, linestyle='-', label='Best-fit model (IC)')
    
    plt.vlines((0.05*u.deg*cluster.cosmo.kpc_proper_per_arcmin(cluster.redshift)).to_value('kpc'),
               0,1, linestyle=':', color='k', label='CTA PSF (1 TeV)')
    plt.vlines(cluster.R500.to_value('kpc'), 0,1, linestyle='--', color='k', label='$R_{500}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'$\frac{dN}{dSdtd\Omega}$ (cm$^{-2}$ s$^{-1}$ deg$^{-2}$)')
    plt.xlim(10,5e3)
    plt.ylim(1e-15,1e-8)
    plt.legend(fontsize=14)
    plt.savefig(cluster.output_dir+'/Haronic_Gamma_profile.pdf')
