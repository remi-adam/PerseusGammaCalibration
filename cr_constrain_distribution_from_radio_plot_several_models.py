"""
This script is part of the CTA key science project dedicated to Perseus. 
Here, we plot the constraints from the different magnetic field models
all together.
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
    sampler_exist = os.path.exists(cluster.output_dir+'/'+model_case+'_sampler.pkl')
    if sampler_exist:
        with open(cluster.output_dir+'/'+model_case+'_sampler.pkl', 'rb') as f:
            sampler = pickle.load(f)
        print('    Existing sampler: '+cluster.output_dir+'/'+model_case+'_sampler.pkl')

    #---------- MCMC parameters
    ndim = len(par0)
    if nwalkers < 2*ndim:
        print('    nwalkers should be at least twice the number of parameters.')
        print('    nwalkers --> '+str(ndim*2))
        print('')
        nwalkers = ndim*2
            
    #----- Define the MCMC
    if sampler_exist:
        if reset_mcmc:
            print('    Reset MCMC even though sampler already exists')
            sampler.reset()
            pos = mcmc_common.chains_starting_point(par0, 0.1, par_min, par_max, nwalkers)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike,
                                            args=[cluster, radio_data, par_min, par_max, par_gprior,
                                                  fit_index, model_case],
                                            pool=Pool(cpu_count()), moves=moves)
        else:
            print('    Start from already existing sampler')
            pos = sampler.chain[:,-1,:]
    else:
        print('    No pre-existing sampler, start from scratch')
        pos = mcmc_common.chains_starting_point(par0, 0.1, par_min, par_max, nwalkers)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike,
                                        args=[cluster, radio_data, par_min, par_max, par_gprior,
                                              fit_index, model_case],
                                        pool=Pool(cpu_count()), moves=moves)
    
    #----- Run the MCMC
    if run_mcmc:
        print('    Runing '+str(mcmc_nsteps)+' MCMC steps')
        res = sampler.run_mcmc(pos, mcmc_nsteps, progress=True)
    
    #----- Save sampler
    with open(cluster.output_dir+'/'+model_case+'_sampler.pkl', 'wb') as output:
        pickle.dump(sampler, output, pickle.HIGHEST_PROTOCOL)

        
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
# Main function
#========================================

if __name__ == "__main__":

    #========== Parameters
    conf        = 68.0
    Nmc         = 100              # Number of Monte Carlo trials
    fit_index   = False            # Fit the spectral index profile
    app_steady  = True             # Application of steady state losses
    mcmc_burnin = 1500                # number of MCMC burnin points
    basedata    = 'Pedlar1990'     # 'Gitti2002', 'Pedlar1990'
    model_case  = 'Hadronic'       # 'Hadronic' or 'Leptonic'
    mag_case    = ['Taylor2006','Bonafede2010up','Walker2017']
    coll1       = ['darkorange', 'darkgreen', 'darkmagenta']
    coll2       = ['orange', 'green', 'magenta']
    output_dir = []
    output_dir0 = '/sps/cta/llr/radam/PerseusGammaCalib'
    for im in mag_case:    
        output_dir.append(output_dir0+'_'+model_case+'_'+im+'_'+basedata)
    
    #========== Define the cluster model
    cluster = []
    for i in range(len(mag_case)):
        if model_case == 'Hadronic':
            cluster0 = perseus_model_library.default_model(directory=output_dir[i])
            cluster0 = perseus_model_library.set_magnetic_field_model(cluster0, case=mag_case[i])
            cluster0 = perseus_model_library.set_pure_hadronic_model(cluster0, ('density', 1.0), 1e-2, 2.5)
            cluster0.Npt_per_decade_integ = 10
        elif model_case == 'Leptonic':
            cluster0 = perseus_model_library.default_model(directory=output_dir[i])
            cluster0 = perseus_model_library.set_magnetic_field_model(cluster0, case=mag_case[i])
            cluster0 = perseus_model_library.set_pure_leptonic_model(cluster0, ('density', 1.0), 1e-5, 2.0)
            if app_steady: cluster0.cre1_loss_model = 'Steady'
            cluster0.Npt_per_decade_integ = 10
        else:
            raise ValueError('Only Hadronic or Leptonic are possible')
        cluster.append(cluster0)
        
    #========== Data
    print('')
    print('-----> Getting the radio data')
    radio_data = perseus_data_library.get_radio_data(cluster[0].cosmo, cluster[0].redshift, prof_file=basedata)
    
    #========== Starting point
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

    #========== Post analysis
    print('')
    print('-----> Starting the post-analysis')

    #========== Get the chains    
    #----- Restore chains
    sampler = []
    for cl in cluster:
        with open(cl.output_dir+'/'+model_case+'_sampler.pkl', 'rb') as f:
            sampler.append(pickle.load(f))

    #----- Burn in
    param_name2 = ['$X_{CRp} (x10^{-2})$', '$\\eta_{CRp}$', '$\\alpha_{CRp}$', 'Norm']
    param_chains = []
    lnL_chains = []
    dfs = []
    for sam in sampler:
        param_chains.append(sam.chain[:, mcmc_burnin:, :])
        lnL_chains.append(sam.lnprobability[:, mcmc_burnin:])
        list_ch = []
        for ip in range(len(param_name)):
            list_ch.append((sam.chain[:, mcmc_burnin:, ip]).flatten())
        dfs.append(pd.DataFrame(np.array(list_ch).T,columns=param_name2))
    ndim = param_chains[0].shape[2]

    #----- Compare PDFs
    plotting.seaborn_corner(dfs,
                            output_fig=output_dir0+'/'+model_case+'_2D_likelihood.pdf', 
                            ci2d=[0.95, 0.68], ci1d=0.68,
                            truth=None, truth_style='star', labels=mag_case,
                            gridsize=50, linewidth=0.75, alpha=(0.1, 0.3, 1.0), n_levels=30,
                            zoom=0.15, add_grid=True,
                            figsize=(12,10), fontsize=12,
                            cols = [('orange',None,'orange','Oranges'),
                                    ('green',None,'green','Greens'), 
                                    ('magenta',None,'magenta','RdPu')])    
    
    #----- Get the best fit parameters
    param_best = np.zeros((ndim, len(sampler)))
    for j in range(len(sampler)):
        wbest = (lnL_chains[j] == np.amax(lnL_chains[j]))
        for i in range(ndim):
            param_best[i,j] = (((param_chains[j][:,:,i])[wbest])[0])
        
    #----- MC parameters
    param_MC = np.zeros((Nmc, ndim, len(sampler)))
    for j in range(len(sampler)):
        param_flat = param_chains[j].reshape(param_chains[j].shape[0]*param_chains[j].shape[1],
                                             param_chains[j].shape[2])
        Nsample = len(param_flat[:,0])-1
        for i in range(Nmc):
            param_MC[i,:,j] = param_flat[np.random.randint(0, high=Nsample), :] # randomly taken from chains
        
    #========== Data versus model
    prof_best, spec_best, idx_best = [], [], []
    for j in range(len(sampler)):
        # Best-fit
        if model_case == 'Hadronic':
            prof_best0, spec_best0, idx_best0 = model_hadronic(param_best[:,j], cluster[j], radio_data)
        if model_case == 'Leptonic':
            prof_best0, spec_best0, idx_best0 = model_leptonic(param_best[:,j], cluster[j], radio_data)
        prof_best.append(prof_best0)
        spec_best.append(spec_best0)
        idx_best.append(idx_best0)
        
    # MC sampling
    prof_mc = np.zeros((len(prof_best0), Nmc, len(sampler)))
    spec_mc = np.zeros((len(spec_best0), Nmc, len(sampler)))
    if idx_best0 is not np.nan:
        idx_mc  = np.zeros((len(idx_best0), Nmc, len(sampler)))
    else:
        idx_mc  = np.zeros((1, Nmc, len(sampler)))
    for j in range(len(sampler)):
        for imc in range(Nmc):
            if model_case == 'Hadronic':
                prof_mci, spec_mci, idx_mci = model_hadronic(param_MC[imc,:,j], cluster[j], radio_data)
            if model_case == 'Leptonic':
                prof_mci, spec_mci, idx_mci = model_leptonic(param_MC[imc,:,j], cluster[j], radio_data)            
            prof_mc[:,imc,j] = prof_mci
            spec_mc[:,imc,j] = spec_mci
            idx_mc[:,imc,j]  = idx_mci

    # Limits
    prof_u, prof_d, spec_u, spec_d, idx_u, idx_d = [], [], [], [], [], []
    for j in range(len(sampler)):
        prof_u.append(np.percentile(np.array(prof_mc[:,:,j]), 100-(100-conf)/2.0, axis=1)*prof_mci.unit)
        prof_d.append(np.percentile(np.array(prof_mc[:,:,j]), (100-conf)/2.0, axis=1)*prof_mci.unit)
        
        spec_u.append(np.percentile(np.array(spec_mc[:,:,j]), 100-(100-conf)/2.0, axis=1)*spec_mci.unit)
        spec_d.append(np.percentile(np.array(spec_mc[:,:,j]), (100-conf)/2.0, axis=1)*spec_mci.unit)
        
        idx_u.append(np.percentile(np.array(idx_mc[:,:,j]), 100-(100-conf)/2.0, axis=1))
        idx_d.append(np.percentile(np.array(idx_mc[:,:,j]), (100-conf)/2.0, axis=1))

    #----- Spectrum
    fig = plt.figure(0, figsize=(8, 6))
    for j in range(len(sampler)):
        # MC
        for imc in range(Nmc):
            plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_mc[:,imc,j],
                     color=coll1[j], alpha=0.05)
        
        # Limits
        plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_u[j].to_value('Jy'),
                 color=coll1[j], linewidth=2, linestyle='--')
        plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_d[j].to_value('Jy'),
                 color=coll1[j], linewidth=2, linestyle='--')
        plt.fill_between(radio_data['spectrum']['freq'].to_value('MHz'), spec_d[j].to_value('Jy'),
                         spec_u[j], color=coll1[j], alpha=0.2)
        
        # Best fit and data
        plt.plot(radio_data['spectrum']['freq'].to_value('MHz'), spec_best[j],
                 color=coll1[j], linewidth=3, label=mag_case[j])
    
    plt.errorbar(radio_data['spectrum']['freq'].to_value('MHz'), radio_data['spectrum']['flux'].to_value('Jy'),
                 radio_data['spectrum']['error'].to_value('Jy'),
                 marker='o', color='k', linestyle='', label='Data')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('frequency (MHz)')
    plt.ylabel('flux (Jy)')
    plt.legend()
    plt.savefig(output_dir0+'/'+model_case+'_Radio_spectrum.pdf')
    plt.close()

    #----- Profile
    fig = plt.figure(0, figsize=(8, 6))
    for j in range(len(sampler)):
        # MC
        for imc in range(Nmc):
            plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_mc[:,imc,j]*(1*prof_mci.unit).to('Jy arcmin-2'),
                     color=coll1[j], alpha=0.05)
        
        # Limits    
        plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_u[j].to_value('Jy arcmin-2'),
                 color=coll1[j], linewidth=2, linestyle='--')
        plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_d[j].to_value('Jy arcmin-2'),
                 color=coll1[j], linewidth=2, linestyle='--')
        plt.fill_between(radio_data['profile']['radius'].to_value('kpc'), prof_d[j].to_value('Jy arcmin-2'),
                         prof_u[j].to_value('Jy arcmin-2'), color=coll1[j], alpha=0.2)
        
        # Best and data
        plt.plot(radio_data['profile']['radius'].to_value('kpc'), prof_best[j]*(1*prof_mci.unit).to('Jy arcmin-2'),
                 color=coll1[j], linewidth=3, label=mag_case[j])
        
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
    plt.savefig(output_dir0+'/'+model_case+'_Radio_profile.pdf')
    plt.close()
    
    #----- Spectral index
    if fit_index:
        fig = plt.figure(0, figsize=(8, 6))
        for j in range(len(sampler)):
            # MC
            for imc in range(Nmc):
                plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_mc[:,imc,j], color=coll1[j], alpha=0.05)
                
            # Limits
            plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_u[j], color=coll1[j],
                     linewidth=2, linestyle='--')
            plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_d[j], color=coll1[j],
                     linewidth=2, linestyle='--')
            plt.fill_between(radio_data['index']['radius'].to_value('kpc'), idx_d[j], idx_u[j],
                             color=coll1[j], alpha=0.2)
            
            # Best and data    
            plt.plot(radio_data['index']['radius'].to_value('kpc'), idx_best[j], color=coll1[j],
                     linewidth=3, label=mag_case[j])
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
        plt.savefig(output_dir0+'/'+model_case+'_Radio_index.pdf')
        plt.close()

    #========== Implication for gamma rays
    energy = np.logspace(-2,6,100)*u.GeV
    radius = np.logspace(0,4,100)*u.kpc

    #---------- Best-fit
    dN_dEdSdt, dN_dSdtdO, fluxH, dNIC_dEdSdt, dNIC_dSdtdO, fluxIC = [], [], [], [], [], []
    for j in range(len(sampler)):
        if model_case == 'Hadronic':
            cluster0 = perseus_model_library.set_pure_hadronic_model(cluster[j], ('density', param_best[1,j]),
                                                                     param_best[0,j]*1e-2, param_best[2,j])
        if model_case == 'Leptonic':
            cluster0 = perseus_model_library.set_pure_leptonic_model(cluster[j], ('density', param_best[1,j]),
                                                                     param_best[0,j]*1e-5, param_best[2,j])
            if app_steady: cluster0.cre1_loss_model = 'Steady'
     
        # Hadronic
        E, dN_dEdSdt0 = cluster0.get_gamma_spectrum(energy, 
                                                    Rmin=None, Rmax=cluster[0].R500,
                                                    type_integral='cylindrical',
                                                    Rmin_los=None, NR500_los=5.0)
        
        r, dN_dSdtdO0 = cluster0.get_gamma_profile(radius, 
                                                   Emin=50*u.GeV, Emax=100*u.TeV, 
                                                   Energy_density=False, Rmin_los=None, NR500_los=5.0)
        
        fluxH0 = cluster0.get_gamma_flux(Emin=50*u.GeV, Emax=None, Energy_density=False,
                                         Rmax=cluster[0].R500, type_integral='cylindrical').to_value('s-1 cm-2')
        
        # Inverse Compton
        E, dNIC_dEdSdt0 = cluster0.get_ic_spectrum(energy, 
                                                   Rmin=None, Rmax=cluster[0].R500,
                                                   type_integral='cylindrical',
                                                   Rmin_los=None, NR500_los=5.0)
        
        r, dNIC_dSdtdO0 = cluster0.get_ic_profile(radius, 
                                                  Emin=50*u.GeV, Emax=100*u.TeV, 
                                                  Energy_density=False, Rmin_los=None, NR500_los=5.0)
        
        fluxIC0 = cluster0.get_ic_flux(Emin=50*u.GeV, Emax=None, Energy_density=False,
                                       Rmax=cluster[0].R500, type_integral='cylindrical').to_value('s-1 cm-2')

        dN_dEdSdt.append(dN_dEdSdt0)
        dN_dSdtdO.append(dN_dSdtdO0)
        fluxH.append(fluxH0)
        dNIC_dEdSdt.append(dNIC_dEdSdt0)
        dNIC_dSdtdO.append(dNIC_dSdtdO0)
        fluxIC.append(fluxIC0)

    #---------- Monte Carlo sampling: hadronic
    prof_g_mc = np.zeros((len(radius),Nmc,len(sampler)))
    spec_g_mc = np.zeros((len(energy),Nmc,len(sampler)))
    flux_g_mc = np.zeros((1,Nmc,len(sampler)))
    for j in range(len(sampler)):
        for imc in range(Nmc):
            if model_case == 'Hadronic':
                cluster0 = perseus_model_library.set_pure_hadronic_model(cluster[j], ('density',
                                                                                      param_MC[imc,1,j]), 
                                                                         param_MC[imc,0,j]*1e-2, param_MC[imc,2,j])
            if model_case == 'Leptonic':
                cluster0 = perseus_model_library.set_pure_leptonic_model(cluster[j], ('density',
                                                                                      param_MC[imc,1,j]),
                                                                         param_MC[imc,0,j]*1e-5, param_MC[imc,2,j])
                if app_steady: cluster0.cre1_loss_model = 'Steady'
                
            spec_g_mc[:,imc,j] = cluster0.get_gamma_spectrum(energy, Rmin=None, Rmax=cluster[j].R500,
                                                             type_integral='cylindrical',
                                                             Rmin_los=None, NR500_los=5.0)[1]
            
            prof_g_mc[:,imc,j] = cluster0.get_gamma_profile(radius, Emin=50*u.GeV, Emax=100*u.TeV, 
                                                            Energy_density=False, Rmin_los=None, NR500_los=5.0)[1]
            
            flux_g_mc[:,imc,j] = cluster0.get_gamma_flux(Emin=50*u.GeV, Emax=None, Energy_density=False,
                                                         Rmax=cluster[j].R500,
                                                         type_integral='cylindrical').to_value('s-1 cm-2')

    #---------- Limits: hadronic
    dN_dEdSdt_u, dN_dEdSdt_d, dN_dSdtdO_u, dN_dSdtdO_d, fluxH_u, fluxH_d = [], [], [], [], [], []
    for j in range(len(sampler)):
        dN_dEdSdt_u.append(np.percentile(np.array(spec_g_mc[:,:,j]), 100-(100-conf)/2.0, axis=1)*dN_dEdSdt0.unit)
        dN_dEdSdt_d.append(np.percentile(np.array(spec_g_mc[:,:,j]), (100-conf)/2.0, axis=1)*dN_dEdSdt0.unit)
        
        dN_dSdtdO_u.append(np.percentile(np.array(prof_g_mc[:,:,j]), 100-(100-conf)/2.0, axis=1)*dN_dSdtdO0.unit)
        dN_dSdtdO_d.append(np.percentile(np.array(prof_g_mc[:,:,j]), (100-conf)/2.0, axis=1)*dN_dSdtdO0.unit)
        
        fluxH_u.append(np.percentile(np.array(flux_g_mc[:,:,j]), 100-(100-conf)/2.0, axis=1))
        fluxH_d.append(np.percentile(np.array(flux_g_mc[:,:,j]), (100-conf)/2.0, axis=1))

    #---------- Monte Carlo sampling: IC
    prof_ic_mc = np.zeros((len(radius),Nmc,len(sampler)))
    spec_ic_mc = np.zeros((len(energy),Nmc,len(sampler)))
    flux_ic_mc = np.zeros((1,Nmc,len(sampler)))
    for j in range(len(sampler)):
        for imc in range(Nmc):
            if model_case == 'Hadronic':
                cluster0 = perseus_model_library.set_pure_hadronic_model(cluster[j], ('density',
                                                                                      param_MC[imc,1,j]), 
                                                                         param_MC[imc,0,j]*1e-2, param_MC[imc,2,j])
            if model_case == 'Leptonic':
                cluster0 = perseus_model_library.set_pure_leptonic_model(cluster[j], ('density',
                                                                                      param_MC[imc,1,j]),
                                                                         param_MC[imc,0,j]*1e-5, param_MC[imc,2,j])
            if app_steady: cluster0.cre1_loss_model = 'Steady'
            
            spec_ic_mc[:,imc,j] = cluster0.get_ic_spectrum(energy, Rmin=None, Rmax=cluster[j].R500,
                                                           type_integral='cylindrical',
                                                           Rmin_los=None, NR500_los=5.0)[1]
            
            prof_ic_mc[:,imc,j] = cluster0.get_ic_profile(radius, Emin=50*u.GeV, Emax=100*u.TeV, 
                                                          Energy_density=False, Rmin_los=None, NR500_los=5.0)[1]
            
            flux_ic_mc[:,imc,j] = cluster0.get_ic_flux(Emin=50*u.GeV, Emax=None, Energy_density=False,
                                                       Rmax=cluster[j].R500,
                                                       type_integral='cylindrical').to_value('s-1 cm-2')
        
    #---------- Limits: IC
    dNIC_dEdSdt_u, dNIC_dEdSdt_d, dNIC_dSdtdO_u, dNIC_dSdtdO_d, fluxIC_u, fluxIC_d = [], [], [], [], [], []
    for j in range(len(sampler)):
        dNIC_dEdSdt_u.append(np.percentile(np.array(spec_ic_mc[:,:,j]),100-(100-conf)/2.,axis=1)*dNIC_dEdSdt0.unit)
        dNIC_dEdSdt_d.append(np.percentile(np.array(spec_ic_mc[:,:,j]),(100-conf)/2.0,axis=1)*dNIC_dEdSdt0.unit)
        
        dNIC_dSdtdO_u.append(np.percentile(np.array(prof_ic_mc[:,:,j]),100-(100-conf)/2.,axis=1)*dNIC_dSdtdO0.unit)
        dNIC_dSdtdO_d.append(np.percentile(np.array(prof_ic_mc[:,:,j]), (100-conf)/2.0, axis=1)*dNIC_dSdtdO0.unit)
        
        fluxIC_u.append(np.percentile(np.array(flux_ic_mc[:,:,j]), 100-(100-conf)/2.0, axis=1))
        fluxIC_d.append(np.percentile(np.array(flux_ic_mc[:,:,j]), (100-conf)/2.0, axis=1))

    #========== Figure
    #----- Spectrum
    fig = plt.figure(0, figsize=(8, 6))
    for j in range(len(sampler)):
        # MC
        for imc in range(Nmc):
            plt.plot(E.to_value('GeV'), (E**2*spec_g_mc[:,imc,j]*dN_dEdSdt0.unit).to_value('MeV cm-2 s-1'),
                     color=coll1[j], alpha=0.05)
        for imc in range(Nmc):
            plt.plot(E.to_value('GeV'), (E**2*spec_ic_mc[:,imc,j]*dNIC_dEdSdt0.unit).to_value('MeV cm-2 s-1'),
                     color=coll2[j], alpha=0.05)
        
        # Limits
        plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt_u[j]).to_value('MeV cm-2 s-1'), color=coll1[j],
                 linewidth=2, linestyle='--')
        plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt_d[j]).to_value('MeV cm-2 s-1'), color=coll1[j],
                 linewidth=2, linestyle='--')
        plt.fill_between(E.to_value('GeV'), (E**2*dN_dEdSdt_u[j]).to_value('MeV cm-2 s-1'),
                         (E**2*dN_dEdSdt_d[j]).to_value('MeV cm-2 s-1'), color=coll1[j], alpha=0.2)
        
        plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt_u[j]).to_value('MeV cm-2 s-1'),
                 color=coll2[j], linewidth=2, linestyle='--')
        plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt_d[j]).to_value('MeV cm-2 s-1'),
                 color=coll2[j], linewidth=2, linestyle='--')
        plt.fill_between(E.to_value('GeV'), (E**2*dNIC_dEdSdt_u[j]).to_value('MeV cm-2 s-1'),
                         (E**2*dNIC_dEdSdt_d[j]).to_value('MeV cm-2 s-1'), color=coll2[j], alpha=0.2)
    
        # Best fit
        plt.plot(E.to_value('GeV'), (E**2*dN_dEdSdt[j]).to_value('MeV cm-2 s-1'),
                 color=coll1[j], linewidth=3)
        plt.plot(E.to_value('GeV'), (E**2*dNIC_dEdSdt[j]).to_value('MeV cm-2 s-1'),
                 color=coll2[j], linewidth=3, linestyle='-')
        
    plt.fill_between([30, 100e3], [0,0], [1e6,1e6], color='red', alpha=0.1, label='CTA energy range')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (GeV)')
    plt.ylabel(r'$\frac{E^2 dN}{dEdSdt}$ (MeV cm$^{-2}$ s$^{-1}$)')
    plt.xlim(1e-2, 2e5)
    plt.ylim(1e-10, 5e-6)
    plt.legend(fontsize=14)
    plt.savefig(output_dir0+'/'+model_case+'_Gamma_spectrum.pdf')
    plt.close()
    
    #----- Profile
    fig = plt.figure(0, figsize=(8, 6))
    for j in range(len(sampler)):
        # MC
        for imc in range(Nmc):
            plt.plot(r.to_value('kpc'), (prof_g_mc[:,imc,j]*dN_dSdtdO0.unit).to_value('cm-2 s-1 deg-2'),
                     color=coll1[j], alpha=0.05)
        for imc in range(Nmc):
            plt.plot(r.to_value('kpc'), (prof_ic_mc[:,imc,j]*dNIC_dSdtdO0.unit).to_value('cm-2 s-1 deg-2'),
                     color=coll2[j], alpha=0.05)
    
        # Limits
        plt.plot(r.to_value('kpc'), (dN_dSdtdO_u[j]).to_value('cm-2 s-1 deg-2'),
                 color=coll1[j], linewidth=2, linestyle='--')
        plt.plot(r.to_value('kpc'), (dN_dSdtdO_d[j]).to_value('cm-2 s-1 deg-2'),
                 color=coll1[j], linewidth=2, linestyle='--')
        plt.fill_between(r.to_value('kpc'), (dN_dSdtdO_u[j]).to_value('cm-2 s-1 deg-2'),
                         (dN_dSdtdO_d[j]).to_value('cm-2 s-1 deg-2'), color=coll1[j], alpha=0.2)
        
        plt.plot(r.to_value('kpc'), (dNIC_dSdtdO_u[j]).to_value('cm-2 s-1 deg-2'),
                 color=coll2[j], linewidth=2, linestyle='--')
        plt.plot(r.to_value('kpc'), (dNIC_dSdtdO_d[j]).to_value('cm-2 s-1 deg-2'),
                 color=coll2[j], linewidth=2, linestyle='--')
        plt.fill_between(r.to_value('kpc'), (dNIC_dSdtdO_u[j]).to_value('cm-2 s-1 deg-2'),
                         (dNIC_dSdtdO_d[j]).to_value('cm-2 s-1 deg-2'), color=coll2[j], alpha=0.2)
    
        # Best-fit
        plt.plot(r.to_value('kpc'), (dN_dSdtdO[j]).to_value('cm-2 s-1 deg-2'),
                 color=coll1[j], linewidth=3)
        plt.plot(r.to_value('kpc'), (dNIC_dSdtdO[j]).to_value('cm-2 s-1 deg-2'),
                 color=coll2[j], linewidth=3, linestyle='-')
    
    plt.vlines((0.05*u.deg*cluster[0].cosmo.kpc_proper_per_arcmin(cluster[0].redshift)).to_value('kpc'),
               0,1, linestyle=':', color='k', label='CTA PSF (1 TeV)')
    plt.vlines(cluster[0].R500.to_value('kpc'), 0,1, linestyle='--', color='k', label='$R_{500}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'$\frac{dN}{dSdtd\Omega}$ (cm$^{-2}$ s$^{-1}$ deg$^{-2}$)')
    plt.xlim(10,5e3)
    plt.ylim(1e-15,1e-8)
    plt.legend(fontsize=14)
    plt.savefig(output_dir0+'/'+model_case+'_Gamma_profile.pdf')
    plt.close()
    
    #----- Flux
    x1 = np.array(flux_g_mc[0,:,0])*1e12
    x2 = np.array(flux_g_mc[0,:,1])*1e12
    x3 = np.array(flux_g_mc[0,:,2])*1e12
    plotting.seaborn_1d([x1,x2,x3],
                        output_fig=output_dir0+'/'+model_case+'_Gamma_flux_H.pdf',
                        ci=0.68, truth=None,
                        label='$F_{\gamma,\ hadronic}$ ($10^{-12}$ s$^{-1}$ cm$^{-2}$)',
                        gridsize=100, alpha=(0.2, 0.4), 
                        figsize=(10,10), fontsize=12,
                        cols=[(coll1[0], coll1[0], coll1[0]),
                              (coll1[1], coll1[1], coll1[1]),
                              (coll1[2], coll1[2], coll1[2])])
    plt.close("all")

    x1 = np.array(flux_ic_mc[0,:,0])*1e12
    x2 = np.array(flux_ic_mc[0,:,1])*1e12
    x3 = np.array(flux_ic_mc[0,:,2])*1e12
    plotting.seaborn_1d([x1,x2,x3],
                        output_fig=output_dir0+'/'+model_case+'_Gamma_flux_IC.pdf',
                        ci=0.68, truth=None,
                        label='$F_{\gamma,\ IC}$ ($10^{-12}$ s$^{-1}$ cm$^{-2}$)',
                        gridsize=100, alpha=(0.2, 0.4), 
                        figsize=(10,10), fontsize=12,
                        cols=[(coll1[0],coll1[0], coll1[0]),
                              (coll1[1],coll1[1], coll1[1]),
                              (coll1[2],coll1[2], coll1[2])])
    plt.close("all")



