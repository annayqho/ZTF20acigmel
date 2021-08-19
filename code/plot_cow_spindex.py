""" Plot the spectral index of 18cow bands over time """

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Initialize plot
fig,ax = plt.subplots(1,1,figsize=(6,4.5))

# Get data
a = pd.read_table("../data/AT2018cow/radio_lc.dat", delimiter='&')
days = a['Day'].values
tel = a['Tel'].values
fef = a['Flux'].values
freq = a['Freq'].values
f = np.array([val.split('pm')[0] for val in fef]).astype(float)
ef = np.array([val.split('pm')[1] for val in fef]).astype(float)

# Round the days to 0.1
days = np.round(days, 1)

# Only use SMA observations
choose = tel=='SMA'
days = days[choose]
f = f[choose]
freq = freq[choose]
ef = ef[choose]

# Remove limits
days = days[ef<99]
f = f[ef<99]
freq = freq[ef<99]
ef = ef[ef<99]

# Remove negative flux values
days = days[f>0]
freq = freq[f>0]
ef = ef[f>0]
f = f[f>0]

# Remove things that are not 5-sigma detections
det = f/ef>=5
days = days[det]
freq = freq[det]
ef = ef[det]
f = f[det]

# Merge similar frequencies
choose = np.logical_and(freq>=215.5, freq<=218)
freq[choose] = 217
choose = np.logical_and(freq>=230.6, freq<=234.6)
freq[choose] = 232
choose = np.logical_and(freq>=341.5, freq<=349)
freq[choose] = 345

udays = np.unique(days)


def fitfunc(x,a,beta):
    return a*(x**beta)


def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))


def plot(d,alpha,ealpha,nu1,nu2):
    # First few epochs have only these frequencies
    if np.logical_and(nu1==217, nu2==232):
        ax.errorbar(
                d, alpha, yerr=ealpha, fmt='o', c='grey', 
                mfc='white', mec='grey', label='217/232')
    elif np.logical_and(nu1==232, nu2==250.6):
        ax.errorbar(d, alpha, yerr=ealpha, fmt='s', c='orange', label='330.8/346.8')
    elif np.logical_and(nu1==232, nu2==246.6):
        ax.errorbar(
                d, alpha, yerr=ealpha, fmt='s', 
                c='orange', label='330.8/346.8', mec='orange', mfc='white')
    elif np.logical_and(nu1==344.8, nu2==360.8):
        ax.errorbar(d, alpha, yerr=ealpha, fmt='D', c='purple', label='344.8/360.8')
    elif np.logical_and(nu1==341.5, nu2==357.5):
        ax.errorbar(d, alpha, yerr=ealpha, fmt='o', mec='k', c='k', mfc='white', label='341.5/357.5')
    elif np.logical_and(nu1==243.3, nu2==259.3):
        ax.errorbar(d, alpha, yerr=ealpha, fmt='D', c='purple', label='243.3/259.3')
    else:
        print(nu1,nu2)
    #    ax.errorbar(d, alpha, yerr=ealpha, fmt='o', c='lightgrey', label='217/232')


# At every epoch, plot the spectral indices
for d in udays:
    choose = days==d
    x = freq[choose]
    y = f[choose]
    ey = ef[choose]

    if sum(choose)==2:
        alpha = np.log(y[1]/y[0])/np.log(x[1]/x[0])
        ealpha = np.abs((1/np.log(x[1]/x[0])) * (1/(y[0]*y[1])) * (f[0]*ef[1]-f[1]*ef[0]))
        plot(d, alpha, ealpha, x[0], x[1])
    elif sum(choose)>2:
        # Fit a spectral index across the whole band ?
        for ii in np.arange(len(x)-1):
            popt,pcov = curve_fit(fitfunc, x, y, sigma=ey, absolute_sigma=False)
            alpha = popt[1]
            ealpha = np.sqrt(pcov[1,1])
            ax.errorbar(d, alpha, yerr=ealpha, fmt='o', c='k', label='217 to 345')
            
legend_without_duplicate_labels(ax)
ax.set_xlabel("Epoch [d]", fontsize=14)
ax.set_ylabel("Spectral index", fontsize=14)
ax.tick_params(axis='both', labelsize=14)
ax.axhline(y=-1.5, ls='--', c='k')
ax.set_xscale('log')
ax.set_xlim(5,30)
ax.set_ylim(-2.5, 1.5)
ax.set_xticks([6, 10, 20, 30])
ax.set_xticklabels([6, 10, 20, 30])
plt.tight_layout()
#plt.savefig("cow_spindex_time.png", dpi=200)
plt.show()
#plt.close()

