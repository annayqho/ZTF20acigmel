""" Plot the peak flux density versus frequency in the rest-frame """

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from get_radio import *


def fitfunc(x,f0,v0,alpha):
    return f0*(x/v0)**(alpha)


fig,ax = plt.subplots(1,1,figsize=(5.5,4))
islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()
freq_obs[freq_obs==34] = 33
freq_unique = np.unique(freq_obs)

x_all = []
y_all = []
ey_all = []

for freq in freq_unique:
    choose = np.logical_and(freq_obs==freq, islim==False)
    if sum(choose)>0:
        print(freq)
        ind = np.argmax(flux_obs[choose])
        peak_flux = flux_obs[choose][ind]
        peak_eflux = eflux_obs[choose][ind]

        # convert to rest-frame
        x = freq*1.2442
        y = peak_flux/1.2442
        ey = peak_eflux/1.2442
        x_all.append(x)
        y_all.append(y)
        ey_all.append(ey)

        if freq in [94.244, 78.756, 33, 10]:
            ax.errorbar(
                    x, y, ey, fmt='o', mec='k', mfc='k', c='k', ms=8)
        else:
            ax.errorbar(x, y, ey, fmt='o', mec='lightgrey', mfc='white', c='lightgrey')


x_all = np.array(x_all)
y_all = np.array(y_all)
ey_all = np.array(ey_all)
choose = x_all<90

p0 = [1, 60, 1]
popt, pcov = curve_fit(fitfunc, x_all[choose], y_all[choose], sigma=ey_all[choose], absolute_sigma=True, p0=p0)

xplot = np.linspace(7,200)
yplot = fitfunc(xplot,*popt)
plt.plot(xplot, yplot, lw=0.5, ls=':', c='k')

ax.scatter(100,100,marker='o', c='lightgrey', label="Sparse light curve")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Rest Frequency [GHz]", fontsize=16)
plt.ylabel(r"Max Observed $f_\nu$ (mJy)", fontsize=16)
plt.tick_params(axis='both', labelsize=12)
plt.xlim(6, 300)
plt.ylim(0.07, 1)
plt.tight_layout()
plt.legend()
plt.show()
#plt.savefig("peak_flux_with_freq.png", dpi=200)
#plt.close()
