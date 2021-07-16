""" Plot spectral index between various bands as a function of time """

import matplotlib.pyplot as plt
import numpy as np
from get_radio import *

cols = ['k', 'k', '#003f5c', '#2f4b7c', '#665191', '#a05195', '#d45087', '#f95d6a', '#ff7c43', '#ffa600']
cols = cols[::-1]
markers = ['o', '*', 'D', 's', 'D', 'P', '*', 'X', 'o', 's', 'D']
msize = [6, 12, 6, 6, 6, 10, 12, 10, 6, 6, 6]


fig,ax = plt.subplots(1,1,figsize=(6,4.5))
islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()

# merge 33 and 34 GHz
freq_obs[freq_obs==34] = 33

# don't include 18 GHz because it's from the ATels...
discard = np.logical_or(freq_obs==18, flux_obs/eflux_obs<3)
islim = islim[~discard]
tel = tel[~discard]
freq_obs = freq_obs[~discard]
days_obs = days_obs[~discard]
flux_obs = flux_obs[~discard]
eflux_obs = eflux_obs[~discard]

z = 0.2442
z = 0

# Convert to the rest-frame
freq = freq_obs * (1+z)
flux = flux_obs / (1+z)
eflux = eflux_obs / (1+z)
days = days_obs / (1+z)

# For each freq, get the one above it
freq_sorted = np.sort(np.unique(freq))
plot_ind = 0

for ii,nu in enumerate(freq_sorted[:-1]):
    nu_lo = nu
    nu_hi = freq_sorted[ii+1]
    #print("checking freq %s/%s" %(nu_lo/1.2442, nu_hi/1.2442))
    print("checking freq %s/%s" %(nu_lo, nu_hi))

    # arrays to plot
    t = []
    alpha = []
    ealpha = []
    arrow = []

    # all the days when this frequency was observed
    low_crit = freq==nu_lo
    days_lo = days[low_crit]
    days_hi = days[freq==nu_hi]
    fluxes_lo = flux[low_crit]
    efluxes_lo = eflux[low_crit]
    fluxes_hi = flux[freq==nu_hi]
    efluxes_hi = eflux[freq==nu_hi]
    islims_lo = islim[low_crit]
    islims_hi = islim[freq==nu_hi]

    # for each day, check if the high frequency was observed within days/10d
    for jj,day_lo in enumerate(days_lo):
        print("checking day %s" %(day_lo))

        is_coeval = np.abs(days_hi-day_lo)<0.1*day_lo

        if sum(is_coeval)>0:
            print("yes, has a pair")
            
            islim_lo = islims_lo[jj]
            islim_hi = islims_hi[is_coeval][0]

            if np.logical_or(islim_lo==False, islim_hi==False):
                print("at least one is a detection")
                # choose one point
                day_hi = days_hi[is_coeval][0]
                
                # average time
                t.append((day_lo+day_hi)/2)

                # get fluxes and errors
                flux_lo = fluxes_lo[jj]
                flux_hi = fluxes_hi[is_coeval][0]
                eflux_lo = efluxes_lo[jj]
                eflux_hi = efluxes_hi[is_coeval][0]
                #print(flux_lo*1.2442, flux_hi*1.2442, eflux_lo*1.2442, eflux_hi*1.2442)
                print(flux_lo, flux_hi, eflux_lo, eflux_hi)

                index = np.log(flux_lo/flux_hi)/np.log(nu_lo/nu_hi)
                alpha.append(index)
                eindex = np.abs((1/np.log(nu_lo/nu_hi)) * (1/(flux_lo*flux_hi)) * (flux_lo*eflux_hi-flux_hi*eflux_lo))
                ealpha.append(eindex)
                print("measured index is %s +/- %s" %(alpha,ealpha))
    
                if islim_lo==True:
                    arrow.append('up')
                elif islim_hi==True:
                    arrow.append('down')
                else:
                    arrow.append('')

    if np.logical_and(len(t)>0, nu_lo!=130.872):
        t = np.array(t)
        alpha = np.array(alpha)
        ealpha = np.array(ealpha)
        arrow = np.array(arrow)
        print("done generating arrays for this frequency")
        print(t)
        print(alpha)
        print(ealpha)
        print(arrow)
        ax.scatter(
                t, alpha, marker=markers[plot_ind], c=cols[plot_ind], 
                label="%s/%s" %(int(nu_lo), int(nu_hi)))
        ax.plot(t, alpha, c=cols[plot_ind])
        choose = arrow==''
        ax.errorbar(
                t[choose], alpha[choose], yerr=ealpha[choose], ms=msize[plot_ind],
                fmt='%s-' %markers[plot_ind], c=cols[plot_ind], lw=0.5)
        choose = arrow!=''
        for k,tval in enumerate(t[choose]):
            print("plotting a pair with a limit")
            fac = 1
            if arrow[k]=='up':
                print("arrow is up")
                fac = 1
            else:
                print("arrow is down")
                fac = -1
            arrowypos = alpha[choose][k]
            ax.arrow(
                tval, arrowypos, 0, fac*0.5,
                length_includes_head=True, 
                color=cols[plot_ind],
                head_width=tval/15, head_length=0.2)
        plot_ind += 1

ax.set_xscale('log')
#ax.set_ylim(-3.0,2.1)
ax.set_xlabel("Days since 2020 Oct 10.0", fontsize=16)
ax.set_ylabel(r"Spectral index $\alpha$ ($f_\nu \propto \nu^{\alpha}$)", fontsize=16)
plt.legend(
        bbox_to_anchor=(0,1.02,1,1.02),loc='lower left',mode='expand', borderaxespad=0.,
        ncol=5, fontsize=10.5)
ax.tick_params(axis='both', labelsize=14)
ax.set_xticks([20,30,40,60,100])
ax.set_xticklabels([20,30,40,60,100])
plt.tight_layout()
#plt.show()
plt.savefig("spindex_time.png", dpi=200)
plt.close()

