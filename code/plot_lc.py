""" Plot radio/mm light curve """

import matplotlib
from astropy.time import Time
from get_radio import *
from ssa_lc import *

# Font sizes
large = 14
medium = 11
small = 9

# Marker sizes
ms=10

# Redshift for transforming to the rest-frame
z = 0.2442

# Initializing plot
fig,axarr = plt.subplots(3,3,figsize=(8,6), sharex=True, sharey=True)


def plot_all(ax):
    """ Plot all of the frequencies as grey in the background """
    islim, tel, freq, days, flux, eflux = get_data_all()
    plot_freqs = np.unique(freq)
    plot_freqs = plot_freqs[plot_freqs < 140]

    for plot_freq in plot_freqs:
        choose = np.logical_and(freq==plot_freq, islim==False)
        ax.plot(days[choose]/(1+z), flux[choose], c='grey', alpha=0.3)


def plot_panel(ax, choose):
    """ Plot a single panel """
    islim, tel, freq, days, flux, eflux = get_data_all()

    # Plot LC at a single frequency
    plot = np.logical_and(islim==False, choose)
    ax.errorbar(
            days[choose]/(1+z), flux[choose], yerr=eflux[choose],
            fmt='o-', c='k', ms=5)


def plot_tpow(ax, choose, pow, ind=1, tjust='right'):
    """ Plot line of t^2 """
    islim, tel, freq, days, flux, eflux = get_data_all()

    xvals = np.linspace(10,200)
    yvals = flux[choose][ind]*(xvals/days[choose][ind])**pow
    ax.plot(xvals/(1+z),yvals,c='k',lw=0.5,ls='--')
    ax.text(days[choose][ind]/(1+z), flux[choose][ind], 
            '$f_\\nu \propto t^{%s}$' %pow, 
            fontsize=small, ha=tjust, va='bottom')


 
if __name__=="__main__":
    islim, tel, freq, days, flux, eflux = get_data_all()

    for ax in axarr.flatten():
        plot_all(ax)

    ax = axarr[0,0]
    choose = freq==94
    plot_panel(ax, choose)
    ax.text(
            0.95,0.25,'94 GHz',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    ax.text(
            0.95,0.12,'(NOEMA)',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    plot_tpow(ax, choose, 2)
    plot_tpow(ax, choose, -3.5, 5, 'left')

    ax = axarr[0,1]
    choose = freq==79
    plot_panel(ax, choose)
    ax.text(
            0.95,0.25,'79 GHz',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    ax.text(
            0.95,0.12,'(NOEMA)',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    plot_tpow(ax, choose, 2)
    plot_tpow(ax, choose, -3, 5, 'left')

    ax = axarr[0,2]
    choose = freq==36.2
    plot_panel(ax, choose)
    ax.text(
            0.95,0.25,'36 GHz',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    ax.text(
            0.95,0.12,'(VLA)',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    plot_tpow(ax, choose, -2.5, 2, 'left')

    ax = axarr[1,0]
    choose = np.logical_and(freq>26, freq<29)
    plot_panel(ax, choose)
    ax.text(
            0.05,0.12,'26/27 GHz (VLA/ATCA)',fontsize=medium,
            transform=ax.transAxes,ha='left', va='top')
    plot_tpow(ax, choose, 2.5)
    plot_tpow(ax, choose, -4, 6, 'left')

    ax = axarr[1,1]
    choose = freq==17.7
    plot_panel(ax, choose)
    ax.text(
            0.05, 0.12,'17.7 GHz (VLA)',fontsize=medium,
            transform=ax.transAxes,ha='left', va='top')
    plot_tpow(ax, choose, -4, 1, 'left')

    ax = axarr[1,2]
    choose = freq==14.5
    plot_panel(ax, choose)
    ax.text(
            0.05, 0.12,'14.5 GHz (VLA)',fontsize=medium,
            transform=ax.transAxes,ha='left', va='top')

    ax = axarr[2,0]
    choose = freq==12
    plot_panel(ax, choose)
    ax.text(
            0.95, 0.12, '12 GHz (VLA)',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    plot_tpow(ax, choose, 2.5)
    plot_tpow(ax, choose, -2.5, 3, 'left')

    ax = axarr[2,1]
    choose = freq==8
    plot_panel(ax, choose)
    ax.text(
            0.95, 0.12, '8 GHz (VLA)',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    plot_tpow(ax, choose, 1)

    ax = axarr[2,2]
    choose = freq==4
    plot_panel(ax, choose)
    ax.text(
            0.95, 0.12, '4 GHz (VLA)',fontsize=medium,
            transform=ax.transAxes,ha='right', va='top')
    plot_tpow(ax, choose, 1)

    # Formatting
    axarr[0,0].set_yscale('log')
    axarr[0,0].set_ylim(0.015,1.3)
    axarr[0,0].set_xscale('log')
    for ax in axarr[:,0]:
        ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=large)
        ax.tick_params(axis='both', labelsize=large)
    for ax in axarr[2,:]:
        ax.set_xlabel("$\Delta t$ [d]", fontsize=large)
        ax.tick_params(axis='both', labelsize=large)
    #for ax in axarr.flatten():
    #    ax.axvline(x=204)

    # ax.set_xlim(11,150)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_yticks([0.1, 1, 0.05])
    # ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.legend(loc='lower right', fontsize=medium)
    # ax.minorticks_off()
    # 
    # #ax.axvline(x=124/1.2442)
    # 
    # ax = axarr[1]
    # ax.set_ylabel("$10^{-15}$ erg cm$^{-2}$ s$^{-1}$", fontsize=large)
    # ax.set_ylim(2E-1,65)
    # ax.set_yticks([0.2, 1, 10])
    # ax.set_xticks([15, 20, 30, 50, 70])
    # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.legend(loc='upper right', fontsize=medium)
    # ax.minorticks_off()
    # 
    # for ax in axarr:

    # Display
    plt.tight_layout()
    plt.show()
    #plt.savefig("radio_lc.png", dpi=300)
    #plt.close()
