""" Plot spectral index between various bands as a function of time """

import matplotlib.pyplot as plt
import numpy as np
from get_radio import *
from format import *

d = get_format()
cols = d['colors']['11']
markers = ['o', '*', 'D', 's', 'D', 'P', '*', 'X', 'o', 's', 'D']
msize = [6, 12, 6, 6, 6, 10, 12, 10, 6, 6, 6]

fig,ax = plt.subplots(1,1,figsize=(3.5,4))
islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()

# merge 33 and 34 GHz
freq_obs[freq_obs==34] = 33

# merge 227 and 230 GHz
freq_obs[freq_obs==226.744] = 230

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

for ii,nu in enumerate(freq_sorted[0:-1]):
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

    # for each day, 
    for jj,day_lo in enumerate(days_lo):
        print("checking day %s" %(day_lo))

        # check if the higher frequency in the pair was observed within dt/10 days
        is_coeval = np.abs(days_hi-day_lo)<0.10*day_lo

        # if so, then...
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

    # Exclude 131 GHz becauase it's off the chart and has a large error bar
    if np.logical_and(len(t)>0, nu_lo!=130.872):
        t = np.array(t)
        alpha = np.array(alpha)
        ealpha = np.array(ealpha)
        arrow = np.array(arrow)
        print("done generating arrays for this frequency")

        # Plot the detections
        choose = arrow==''
        ax.errorbar(
                t[choose], alpha[choose], yerr=ealpha[choose], ms=msize[plot_ind],
                fmt='%s-' %markers[plot_ind], c=cols[plot_ind], lw=2,
                elinewidth=0.5)

        # Plot something for the legend
        ax.scatter(
                -100, -100, marker=markers[plot_ind], c=cols[plot_ind],
                label="%s/%s" %(int(nu_lo), int(nu_hi)))

        # Plot the non-detections
        choose = arrow!=''
        ax.errorbar(
                t[choose], alpha[choose], yerr=[0]*len(t[choose]), ms=msize[plot_ind],
                fmt='%s' %markers[plot_ind], c=cols[plot_ind], lw=2)

        # Connect the non-detections to the detections
        if nu_lo==10:
            ax.plot(t[0:3], alpha[0:3], c=cols[plot_ind], lw=2, ls='--')
        elif nu_lo==22:
            ax.plot(t[1:], alpha[1:], c=cols[plot_ind], lw=2, ls='--')
        elif nu_lo==78.756:
            ax.plot(t[4:], alpha[4:], c=cols[plot_ind], lw=2, ls='--')

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

ax.axhline(y=-1.5, lw=1, c='grey', ls=':')

ax.set_xscale('log')
ax.set_ylim(-3,3)
ax.set_xlabel("Days since 2020 Oct 10.0", fontsize=d['font_med'])
ax.set_ylabel(r"Spectral index $\beta$ ($f_\nu \propto \nu^{\beta}$)", 
        fontsize=d['font_med'])
plt.legend(
        bbox_to_anchor=(0,1.02,1,1.02),loc='lower left',mode='expand', 
        borderaxespad=0., ncol=3, fontsize=d['font_small'])
#plt.legend(
#        bbox_to_anchor=(0,1.02,1,1.02),loc='lower left',mode='expand', 
#        ncol=4, fontsize=d['font_small'])
ax.tick_params(axis='both', labelsize=d['font_med'])
ax.set_xticks([20,30,40,60,100])
ax.set_xticklabels([20,30,40,60,100])
plt.tight_layout()
plt.show()
#plt.savefig("spindex_time.png", dpi=200, bbox_inches='tight')
#plt.close()

