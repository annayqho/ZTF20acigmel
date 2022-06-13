#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 13:46:50 2021

@author: yuhanyao
"""
import numpy as np

import corner
import emcee

from matplotlib.ticker import MaxNLocator

import matplotlib
import matplotlib.pyplot as plt
fs= 10
matplotlib.rcParams['font.size']=fs
ms = 6
matplotlib.rcParams['lines.markersize']=ms


all_pars = ["lgL_peak", "nu_peak", "beta1", "beta2", "lgs", "lgsigma0"]
prior_dict0 = {par:{} for par in all_pars}


def double_PL_smooth(p, nu):
    """
    nu is the frequency in transient's rest-frame
    """
    lgL_peak = p[0]   # Lnu at peak
    L_peak = 10**lgL_peak
    nu_peak = p[1]
    beta1 = p[2]
    beta2 = p[3]
    lgs = p[4]
    s = 10**lgs
    
    Lnu_model = L_peak * ((nu/nu_peak)**(-1*s*beta1) + (nu/nu_peak)**(-1*s*beta2))**(-1/s)
    nuLnu_model = Lnu_model * nu
    return nuLnu_model # erg/s
    

def get_default_priordict(nu, nuLnu, nuLnu_unc):
    prior_dict_defaults = {}
    
    prior_dict_defaults["lgL_peak"] =  {"min": np.log10(max(nuLnu/nu))-1, 
                                        "max": np.log10(max(nuLnu/nu))+1, 
                                        "sigma": None, 
                                        "value": np.log10(max(nuLnu/nu))
                                        }
    
    prior_dict_defaults["nu_peak"] = {"min": 3,
                                      "max": 100,
                                      "sigma": None,
                                      "value": 10
                                      }
    
    prior_dict_defaults["beta1"] = {"min": 0,
                                    "max": 5,
                                    "sigma": None,
                                    "value": 2.5
                                    }
    
    prior_dict_defaults["beta2"] = {"min": -5,
                                    "max": 0,
                                    "sigma": None,
                                    "value": -1
                                    }
    
    prior_dict_defaults["lgs"] = {"min": -2,
                                "max": 1,
                                "sigma": None,
                                "value": 0
                                }
    
    prior_dict_defaults["lgsigma0"] = {"min": np.log10(max(nuLnu))-5,
                                       "max": np.log10(max(nuLnu))-0.5,
                                       "sigma": None,
                                       "value": np.log10(max(nuLnu))-3
                                       }
    
    return prior_dict_defaults
    
    
    



def run_radio_fit(y, yerr, nu, prior_dict, 
                  nsteps = 1000, burnin = 500,
                  plotname = None):
    """
    y: nuLnu
    yerr: uncertainty of nuLnu
    nu is the frequency in transient's rest-frame
    """
    nwalkers = 100
    
    print ('# of datapoints used in MCMC fitting: %d'%len(nu))
    
    prior_dict_defaults = get_default_priordict(nu, y, yerr)
    
    par_names = np.array(["lgL_peak", "nu_peak", "beta1", "beta2", "lgs", "lgsigma0"])
    def model_func(p, nu):
        return double_PL_smooth(p, nu)
    
    # overwrite input parameters with defaults
    for par in par_names:
        for key in ("min", "max", "value", "sigma"):
            if prior_dict[par].get(key) is None:
                prior_dict[par][key] = prior_dict_defaults[par][key]
                print ("%6s: "%par, "using default value for %6s: %s"%(key, prior_dict[par][key]))
           
    for par in par_names:
        if prior_dict[par].get("sigma") is not None:
            prior_dict[par]["norm"] = -np.log(prior_dict[par]["sigma"])
            prior_dict[par]["ivar"] = 1/prior_dict[par]["sigma"]**2
            
    print ("full prior dict:")
    for par in par_names:
        print ("  ", par, prior_dict[par])
        
        
    # Set up the sampler.
    ndim = len(par_names)
    guess_pos = [prior_dict[par]["value"] for par in par_names]
    
    for i in range(len(par_names)):
        par = par_names[i]
        if (guess_pos[i]<prior_dict[par]["min"]) or (guess_pos[i]>prior_dict[par]["max"]):
            print ('WARNING!! the start position ({0:0.3f}) is outside the limits of this prior:'.format(guess_pos[i]), prior_dict[par])
        
    # Prior function 
    par_with_gauss_prior = [par for par in par_names if prior_dict[par]['sigma']]
    def lnprior(theta):
        # check limits
        for i in range(len(par_names)):
            par = par_names[i]
            if theta[i]<prior_dict[par]["min"]:
                return -np.inf
            if theta[i]>prior_dict[par]["max"]:
                return -np.inf    
        # apply Gausian priors (if any)
        out = 0
        for i in range(len(par_names)):
            par = par_names[i]
            if par in par_with_gauss_prior:
                term1 = -0.5*(theta[i]-prior_dict[par]["value"])**2 * prior_dict[par]["ivar"]
                out += term1 + prior_dict[par]["norm"]
        return out


    # Define the likelihood
    def lnlike(theta, y, yerr, nu):
        # get whatever the model is
        model_y = model_func(theta[0:-1], nu)
        
        lgsigma0 = theta[-1]
        sigma0 = 10**lgsigma0
        sigma2 = yerr**2 + sigma0**2
        
        chi2_term = -0.5 * np.sum( (y-model_y)**2/sigma2)
        error_term = np.sum(np.log(1/np.sqrt(2*np.pi*sigma2)))
        lnL = chi2_term + error_term
        return lnL

    # Define the probability function as likelihood * prior.
    def lnprob(theta, y, yerr, nu):        
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        out = lp + lnlike(theta, y, yerr, nu)
        return out


    # Set up the sampler.
    print("MCMC start:")
    for l in range(ndim):
        print ("{0:7}  = {1:0.2f}".format(par_names[l], guess_pos[l])) 
    
    # add small amount of scatter to start of walkers
    pos = [guess_pos + 1e-3*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(y, yerr, nu))

    # Clear and run the production chain.    
    sampler.run_mcmc(pos, nsteps, rstate0=np.random.get_state(), progress=True)
    idx_par_pick = np.ones(len(par_names), dtype = bool)
    
    for par in prior_dict:
        if "sigma" in prior_dict[par]:
            if prior_dict[par]["sigma"] == 1e-5:
                print ("  %s is fixed --> skip plotting... "%par)
                ind = par_names==par
                idx_par_pick[ind] = False
    
    par_names_walk = par_names[idx_par_pick]
            
    if plotname:  
        fig, axes = plt.subplots(len(par_names_walk), 1, sharex=True, figsize=(10.4, 11.7))

        for l, par in enumerate(par_names_walk):
            idx = np.where(par_names==par)[0][0]
            axes[l].plot(sampler.chain[:, :, idx].T, color="k", alpha=0.4)
            axes[l].yaxis.set_major_locator(MaxNLocator(5))
            axes[l].set_ylabel(par)

        axes[l].set_xlabel("step number")

        fig.tight_layout(h_pad=0.0)
        fig.savefig(plotname+"-mcmc-walkers.png")
        plt.close()

    # Make the corner plot
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

    if plotname:    
        fig = corner.corner(samples[:,idx_par_pick], labels=par_names_walk)
        fig.savefig(plotname+"-mcmc-triangle.png")
        plt.close()

    # Compute the quantiles.
    credit_int = 68.2689492137086 # used for mcmc intervals
    out_tup = [np.array([v[1], v[1]-v[0], v[2]-v[1]]) for v in zip(*np.percentile(samples, [50-credit_int/2., 50, 50+credit_int/2.], axis=0))] #[16, 50, 84]
    out_dict = {par:out_tup[l] for l, par in enumerate(par_names)}

    # setup defaults
    xl = np.arange(0.1, 35, 0.1) # x-axis of output radio SED
    
    # Plot some samples onto the data.    
    Nlc = min(200, len(samples[:,0]))
    
    # we compute both the luminosity at nu_kc, the BB luminosity, radius, and temperature
    mcmc_lcs = np.zeros((Nlc, len(xl)))
    
    selected_pars = samples[np.random.randint(len(samples), size=Nlc)]
    for l in range(Nlc):
        parms = selected_pars[l]       
        
        lc = model_func(parms[0:-1], xl)

        # move this into the else-statement if we want to reject walkers
        mcmc_lcs[l, :] = lc 
        
        # parms[1]: lg(Lpeak) of this walker
        if  abs(parms[1]-out_tup[1][0])>5*out_tup[1][1] or \
            abs(parms[-3]-out_tup[-3][0])>5*out_tup[-3][1]:
            
            print('troublesome sample?: {0[0]:0.2f}, {0[1]:0.2f}, {0[2]:0.2f}, {0[3]:0.2f}, {0[4]:0.2f} {0[5]:0.2f}'.format(parms[0:6]))
    

    print("MCMC result:")
    for par in par_names:
        print ("{0:10} = {1[0]:0.3f} -{1[1]:0.3f} +{1[2]:0.3f}".format(par, out_dict[par])) 
    print ("")

    out_dict["lc"] = {'nu':xl, 
                      'nuLnu':mcmc_lcs}
    
    return out_dict 