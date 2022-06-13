""" Script from Yuhan Yao, modified by Anna Ho """

import os
import pickle
import numpy as np
import astropy.constants as const
import sys
sys.path.append("/Users/annaho/Dropbox/astronomy/papers_active/ZTF20acigmel/code")
from get_radio import *
from basics import get_z


from astropy.table import Table
import astropy.io.ascii as asci

from astropy.cosmology import Planck15
import matplotlib
import matplotlib.pyplot as plt
fs= 10
matplotlib.rcParams['font.size']=fs

from pl_models import prior_dict0, run_radio_fit

z = get_z()
D = Planck15.luminosity_distance(z).cgs.value


def plot_models(ax, mcmc_mod1):
    color = "k"
    lc = mcmc_mod1["lc"]
    nu = lc["nu"]
    nuLnus = lc["nuLnu"]
    for l in range(len(nuLnus)):
        # for 100 draws of posterior
        if l%2 ==0:
            nuLnu = nuLnus[l] # bolometric luminosity as a function of time
            ax.plot(nu, nuLnu / nu, '-', color=color, 
                    alpha=0.03, zorder=1)
    
    
def custom_ax(ax):
    ax.semilogx()
    ax.semilogy()
    xmin = 1
    xmax = 35
    ax.set_xlim(xmin, xmax)
    ymin = 1e+28
    ymax = 3e+29
    ax.set_ylim(ymin, ymax)
    
    
def mcmc_model_fit_e1():
    filename = "../../data/data_20mrf/radio/radio_e1.dat"
    tb = asci.read(filename)
    
    nu_obs = tb["nu"].data
    nu = nu_obs * (1+z)
    
    fnu_obs = tb["fnu"].data * 1e-6 * 1e-23
    fnu_obs_unc = tb["fnu_unc"].data * 1e-6 * 1e-23
    
    nuLnu = (4 * np.pi * D_cm**2) * fnu_obs * nu_obs
    nuLnu_unc = (4 * np.pi * D_cm**2) * fnu_obs_unc * nu_obs
    
    #plt.errorbar(nu_obs, nuLnu/nu, nuLnu_unc/nu, fmt = ".k")
    #plt.semilogx()
    #plt.semilogy()
    
    prior_dict = prior_dict0.copy()
    for k in prior_dict:
        for j in prior_dict[k]:
            prior_dict[k][j] = None
            
    prior_dict["lgL_peak"]["value"] = 29.17
    prior_dict["lgL_peak"]["min"] = 28.5
    prior_dict["lgL_peak"]["max"] = 29.5
    
    prior_dict["nu_peak"]["value"] = 7
    prior_dict["nu_peak"]["min"] = 5
    prior_dict["nu_peak"]["max"] = 9
    
    prior_dict["beta1"]["value"] = 2.5
    prior_dict["beta1"]["min"] = 1
    prior_dict["beta1"]["max"] = 5
    
    prior_dict["beta2"]["value"] = -1
    prior_dict["beta2"]["min"] = -2
    prior_dict["beta2"]["max"] = 0
    
    mcmc_dir = "./mcmc_result"
    name = "AT2020mrf"
    
    print ("")
    print ("Model 1")
    print ("==============================================================================")
    print ("Fix beta1 at 2.5 (optically thick spectral index) ")
    print ("and beta2 (optically thin spectral index) at -1")
    print ("==============================================================================")
    
    prior_dict["beta1"]["sigma"] = 1e-5
    prior_dict["beta2"]["sigma"] = 1e-5
    
    prior_dict["lgsigma0"]["value"] = 1
    prior_dict["lgsigma0"]["max"] = 2
    prior_dict["lgsigma0"]["min"] = 0
    prior_dict["lgsigma0"]["sigma"] = 1e-5
    
    fit_name = 'e1_model1'
    fname = '%s/%s_%s.pickle'%(mcmc_dir, name, fit_name)
    
    if not os.path.isfile(fname):
        plotname = fname.replace('.pickle', '')
        nsteps = 2000
        mcmc_mod1 = run_radio_fit(nuLnu, nuLnu_unc, nu, prior_dict, 
                                  nsteps = nsteps, burnin = nsteps-500,
                                  plotname = plotname)
        pickle.dump( (mcmc_mod1, prior_dict), open(fname,'wb'))
    else:
        print ('reading', fname)
        mcmc_mod1 = pickle.load(open(fname,'rb'),encoding='latin1')[0]
        
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    ax.errorbar(nu, nuLnu/nu, nuLnu_unc/nu, fmt = ".k")
    plot_models(ax, mcmc_mod1)
    custom_ax(ax)
    plt.tight_layout()
    plt.savefig(fname.replace('.pickle', '_lc.png'))
    plt.close()
    
    print ("")
    print ("Model 2")
    print ("==============================================================================")
    print ("Fix beta1 at 2.5 (optically thick spectral index) ")
    print ("and free beta2 (optically thin spectral index)")
    print ("==============================================================================")
    
    prior_dict["beta2"]["sigma"] = None
    
    prior_dict["lgsigma0"]["value"] = 1
    prior_dict["lgsigma0"]["max"] = 2
    prior_dict["lgsigma0"]["min"] = 0
    prior_dict["lgsigma0"]["sigma"] = 1e-5
    
    fit_name = 'e1_model2'
    fname = '%s/%s_%s.pickle'%(mcmc_dir, name, fit_name)
    
    if not os.path.isfile(fname):
        plotname = fname.replace('.pickle', '')
        nsteps = 2000
        mcmc_mod2 = run_radio_fit(nuLnu, nuLnu_unc, nu, prior_dict, 
                                  nsteps = nsteps, burnin = nsteps-500,
                                  plotname = plotname)
        pickle.dump( (mcmc_mod2, prior_dict), open(fname,'wb'))
    else:
        print ('reading', fname)
        mcmc_mod2 = pickle.load(open(fname,'rb'),encoding='latin1')[0]
        
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    ax.errorbar(nu, nuLnu/nu, nuLnu_unc/nu, fmt = ".k")
    plot_models(ax, mcmc_mod2)
    custom_ax(ax)
    plt.tight_layout()
    plt.savefig(fname.replace('.pickle', '_lc.png'))
    plt.close()
    
    print ("")
    print ("Model 3")
    print ("==============================================================================")
    print ("Free beta1 (optically thick spectral index) ")
    print ("and free beta2 (optically thin spectral index)")
    print ("==============================================================================")
    
    prior_dict["beta1"]["sigma"] = None
    prior_dict["beta2"]["sigma"] = None
    
    prior_dict["lgsigma0"]["value"] = 1
    prior_dict["lgsigma0"]["max"] = 2
    prior_dict["lgsigma0"]["min"] = 0
    prior_dict["lgsigma0"]["sigma"] = 1e-5
    
    fit_name = 'e1_model3'
    fname = '%s/%s_%s.pickle'%(mcmc_dir, name, fit_name)
    
    if not os.path.isfile(fname):
        plotname = fname.replace('.pickle', '')
        nsteps = 2000
        mcmc_mod3 = run_radio_fit(nuLnu, nuLnu_unc, nu, prior_dict, 
                                  nsteps = nsteps, burnin = nsteps-500,
                                  plotname = plotname)
        pickle.dump( (mcmc_mod3, prior_dict), open(fname,'wb'))
    else:
        print ('reading', fname)
        mcmc_mod3 = pickle.load(open(fname,'rb'),encoding='latin1')[0]
        
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    ax.errorbar(nu, nuLnu/nu, nuLnu_unc/nu, fmt = ".k")
    plot_models(ax, mcmc_mod3)
    custom_ax(ax)
    plt.tight_layout()
    plt.savefig(fname.replace('.pickle', '_lc.png'))
    plt.close()
    

if __name__=="__main__":
    islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()

    # Put into the rest-frame
    freq = freq_obs[islim==False] * (1+z)
    flux = flux_obs[islim==False] / (1+z)
    eflux = eflux_obs[islim==False] / (1+z)
    days = days_obs[islim==False] / (1+z)

    # 

    mcmc_model_fit_e1()
    
    
    
    




