""" Plot the SED of AT2020xnd on a few different days """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from get_radio import *
from format import *


if __name__=="__main__":
    d = get_format()
    cols = d['colors']['9']
    cols = cols[::-1]
    markers = ['o', 's', 'D', 'P', '*', 'X', 'o', 's', 'D']
    msize = [6, 6, 6, 10, 12, 10, 6, 6, 6]

    fig,ax = plt.subplots(1, 1, figsize=(5.0,3.5))
    islim, tel, freq, days, flux, eflux = get_data_all()

    bins = [18, 24, 30.3, 38, 46, 51.9, 71, 95, 132]

    # For each time bin,
    for b,bin in enumerate(bins):
        col = cols[b]

        # choose all points within a dt/10d window
        choose = np.logical_and.reduce((
            days>bin-bin/20, days<bin+bin/20))#, eflux>0))

        # sort in order of frequency
        order = np.argsort(freq[choose])

        # plot the detections
        det = eflux[choose][order]>0
        ax.errorbar(
                freq[choose][order][det], flux[choose][order][det], 
                eflux[choose][order][det], 
                fmt='%s' %markers[b], c=col, label=None, ms=msize[b], lw=2,
                elinewidth=0.5)
        # plot the non-detections
        ax.errorbar(
                freq[choose][order][~det], flux[choose][order][~det], 
                eflux[choose][order][~det], fmt='%s' %markers[b], 
                mec=col, mfc='white', label=None, ms=msize[b], lw=2,
                elinewidth=0.5)
        ax.scatter(0, 0, marker='%s' %markers[b], c=col, label='%s d' %bin)

    # Formatting
    ax.set_ylabel("Flux Density [mJy]", fontsize=d['font_med'])
    ax.set_xlabel("Observed Frequency [GHz]", fontsize=d['font_med'])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.03, 1.2)
    ax.set_xticks([5,10,20,50,100,200])
    ax.set_xticklabels([5,10,20,50,100,200])
    ax.set_yticks([0.05, 0.1, 0.2, 0.5, 1])
    ax.set_yticklabels([0.05, 0.1, 0.2, 0.5, 1])
    ax.set_xlim(5, 230)
    ax.tick_params(axis='both', labelsize=d['font_small'])
    ax.legend(loc='lower right', fontsize=d['font_small'], ncol=2)

    # Final formatting
    plt.tight_layout()
    plt.show()
    #plt.savefig("radio_sed.png", dpi=300)
    #plt.close()



