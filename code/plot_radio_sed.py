""" Plot the SED of AT2020xnd on a few different days """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from get_radio import *


if __name__=="__main__":
    cols = ['blue', 'green', 'yellow', '#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600', 'lightgrey']

    # for 9 points
    cols = ['k', '#003f5c', '#2f4b7c', '#665191', '#a05195', '#d45087', '#f95d6a', '#ff7c43', '#ffa600']
    cols = cols[::-1]
    markers = ['o', 's', 'D', 'P', '*', 'X', 'o', 's', 'D']
    msize = [6, 6, 6, 10, 12, 10, 6, 6, 6]

    fig,ax = plt.subplots(1, 1, figsize=(6,4.5))
    islim, tel, freq, days, flux, eflux = get_data_all()

    bins = [18, 24, 30.3, 38, 46, 51.9, 71, 95, 132]
    #bins = bins[6:7]

    for b,bin in enumerate(bins):
        col = cols[b]
        choose = np.logical_and.reduce((
            days>bin-bin/20, days<bin+bin/20, eflux>0))
        #if bin!=71:
        if bin!=100:
            order = np.argsort(freq[choose])
            ax.errorbar(
                    freq[choose][order], flux[choose][order], eflux[choose][order], 
                    fmt='%s-' %markers[b], c=col, label=None, ms=msize[b], lw=2,
                    elinewidth=0.5)
            ax.scatter(0, 0, marker='%s' %markers[b], c=col, label='%s d' %bin)
            nondet = np.logical_and(choose, eflux==0)
            #for i,x in enumerate(freq[nondet]):
            #    print(x)
            #    y = flux[nondet][i]
            #    print(y)
            #    ax.arrow(x, y, 
            #            0, -y/2, head_width=x/7, 
            #            length_includes_head=True, head_length=y/7, color=col)

        else:
            #keep = freq[choose] < 80
            keep = freq[choose] < 500
            ax.errorbar(
                    freq[choose][keep], flux[choose][keep], eflux[choose][keep], 
                    fmt='%s-' %markers[b], c=col, label=None, ms=msize[b], lw=2,
                    elinewidth=0.5)
            ax.scatter(0, 0, marker='%s' %markers[b], c=col, label='%s d' %bin)


    #day17(days, freq, flux, eflux, ax, cols[-1])
    #day24(axarr[0])
    #day36(axarr[1])
    #day46(days, freq, flux, eflux, axarr[2])
    #day71(days, freq, flux, eflux, axarr[3])
    #day94(days, freq, flux, eflux, axarr[3,0])
    #day131(days, freq, flux, eflux, axarr[3,1])

    # Formatting
    ax.set_ylabel("Flux Density [mJy]", fontsize=14)
    ax.set_xlabel("Observed Frequency [GHz]", fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.03, 1.2)
    ax.set_xlim(5, 230)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(loc='lower right', fontsize=11, ncol=2)

    # Final formatting
    plt.tight_layout()
    #plt.show()
    plt.savefig("radio_sed.png", dpi=300)
    plt.close()



