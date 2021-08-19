""" Plot the SED of AT2020xnd on a few different days """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from get_radio import *


def interpolate(target_day, x1, y1, ey1):
    """ general function to interpolate """
    print("I will interpolate this flux array")
    print(x1)
    print(y1)
    print("on this day")
    print(target_day)
    nsim = 1000
    newflux = np.zeros(nsim)
    ysamples = np.zeros((nsim,len(y1)))
    for ii,val in enumerate(y1):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey1[ii],size=nsim)
    for ii in np.arange(nsim):
        ysample = ysamples[ii]
        newflux[ii] = np.interp(target_day, x1, ysample)
    y = np.mean(newflux)
    ey = np.std(newflux)
    print(y, ey)
    return y, ey


def day17(days, freq, flux, eflux, ax):
    """
    Plot the VLA and NOEMA observations from Day 17
    """

    day = 17
    # Plot the VLA LC
    x = 8
    choose = np.logical_and(days==day,freq==x)
    y = np.min(flux[choose])
    ax.scatter(x, y, marker='_', c='k')
    ax.arrow(x, y, 0, -y/4, head_width=x/10, 
            length_includes_head=True, head_length=y/10, color='k')
    
    x = 12 
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    # Plot the NOEMA LC
    x = 79
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 94
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=17d", transform=ax.transAxes)


def day24(ax):
    """ 
    Plot VLA and NOEMA and SMA upper limit
    """

    choose = np.logical_and(days>20, days<26)
    order = np.argsort(freq[choose])

    # VLA
    ax.errorbar(
            freq[choose][order][0:3], flux[choose][order][0:3], eflux[choose][order][0:3], 
            fmt='o', c='k', label="VLA (25d)")

    # NOEMA
    ax.errorbar(
            freq[choose][order][3:5], flux[choose][order][3:5], eflux[choose][order][3:5], 
            fmt='s', c='k', label="NOEMA (24d)")
   
    # SMA
    x = freq[choose][order][5]
    y = flux[choose][order][5]
    ax.scatter(x, y, 
            marker='D', c='k', label="SMA (21d)")
    ax.arrow(x, y, 
            0, -y/2, head_width=x/7, 
            length_includes_head=True, head_length=y/7, color='k')

    # repeat for 17d
    choose = np.logical_and(days==17, freq==12)

    # VLA
    ax.errorbar(
            freq[choose], flux[choose], eflux[choose], 
            fmt='o', mec='k', mfc='white', label="VLA (17d)", c='k')

    choose = np.logical_and(days==17, freq>70)

    # NOEMA
    ax.errorbar(
            freq[choose], flux[choose], eflux[choose], 
            fmt='s', c='k', mec='k', mfc='white', label="NOEMA (17d)")

    ax.legend(loc='lower right', ncol=2, fontsize=8)

    # Plot the nu^{-0.5} line
    # xplot = np.linspace(30,200)
    # yplot = y*(xplot/x)**(-0.55)
    # ax.plot(xplot, yplot, c='k', lw=0.5, ls='--')
    # ax.text(
    #         0.99, 0.99, r'$F_\nu \propto \nu^{-0.55\pm0.15}$', 
    #         ha='right', va='top', transform=ax.transAxes)


    # # Plot the nu^2.5 line
    # xplot = np.linspace(10,100)
    # yplot = y*(xplot/x)**2.5
    # ax.plot(xplot, yplot, c='k', lw=0.5, ls='--')
    # ax.text(x, y, r'$F_\nu \propto \nu^{5/2}$', ha='right')

    # Label
    # ax.text(0.05,0.85,"$\Delta t$=21-25d", transform=ax.transAxes)


def day31(ax):
    """ 
    Plot NOEMA from Day 31, ATCA from Day 28
    """
    choose = np.logical_and(days>27, days<32)
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='o', c='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=28-31d", transform=ax.transAxes)


def day36(ax):
    """ 
    Plot VLA from day 35 and NOEMA from day 39
    """
    choose = np.logical_and(days>34, days<40)
    order = np.argsort(freq[choose])

    x = freq[choose][order][0:3]
    y = flux[choose][order][0:3]
    ey = eflux[choose][order][0:3]
    ax.errorbar(x, y, ey, fmt='o', c='k', label="VLA (36d)")

    x = freq[choose][order][3:5]
    y = flux[choose][order][3:5]
    ey = eflux[choose][order][3:5]
    ax.errorbar(x, y, ey, fmt='s', c='k', label="NOEMA (39d)")

    x = freq[choose][order][5]
    y = flux[choose][order][5]
    ey = eflux[choose][order][5]
    ax.errorbar(x, y, ey, fmt='D', c='k', label="SMA (35d)")
    ax.arrow(x, y, 
            0, -y/2, head_width=x/7, 
            length_includes_head=True, head_length=y/7, color='k')

    # Label
    #ax.text(0.05,0.85,"$\Delta t$=36-39d", transform=ax.transAxes)
    ax.legend(loc='lower right', fontsize=8)


def day46(days, freq, flux, eflux, ax):
    """ 
    Interpolate VLA/ATCA, fix NOEMA
    """
    choose = days==51
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='o', c='k',
            label="VLA (51d)")

    choose = days==46
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='s', c='k',
            label="NOEMA (46d)")

    # SMA upper limit
    choose = days==47
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='D', c='k',
            label="SMA (47d)")
    ax.arrow(freq[choose][0], flux[choose][0], 
            0, -flux[choose][0]/2, head_width=freq[choose][0]/6, 
            length_includes_head=True, head_length=flux[choose][0]/5,
            color='k')

    # Label
    #ax.text(0.05,0.85,"$\Delta t$=46-51d", transform=ax.transAxes)
    ax.legend(ncol=2, fontsize=8)


def day71(days, freq, flux, eflux, ax):
    """ 
    Just plot VLA
    """
    choose = days==71
    ax.errorbar(
            freq[choose], flux[choose], eflux[choose], 
            fmt='o', c='k', label="VLA (71d)")

    choose = days==94
    ax.errorbar(
            freq[choose], flux[choose], eflux[choose], 
            fmt='o', c='grey', label="VLA (94d)")

    choose = days==131
    ax.errorbar(
            freq[choose], flux[choose], eflux[choose], 
            fmt='o', c='k', mec='k', mfc='white', label="VLA (131d)")
    ax.arrow(freq[choose][-1], flux[choose][-1], 
            0, -flux[choose][-1]/2, head_width=freq[choose][-1]/6, 
            length_includes_head=True, head_length=flux[choose][-1]/5,
            color='k')

    choose = days==67
    ax.errorbar(
            freq[choose], flux[choose], eflux[choose], 
            fmt='s', c='k', label="NOEMA (67d)")
    ax.arrow(freq[choose][-1], flux[choose][-1], 
            0, -flux[choose][-1]/2, head_width=freq[choose][-1]/6, 
            length_includes_head=True, head_length=flux[choose][-1]/5,
            color='k')

    # Label
    #ax.text(0.05,0.85,"$\Delta t$=67-71d", transform=ax.transAxes)
    ax.legend(fontsize=8, loc='upper right', ncol=4, columnspacing=0.5)


def day94(days, freq, flux, eflux, ax):
    """ 
    Just plot VLA
    """
    choose = days==94
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='o', c='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=94d", transform=ax.transAxes)


def day131(days, freq, flux, eflux, ax):
    """ 
    Just plot VLA
    """
    choose = days==131
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='o', c='k')
    ax.arrow(freq[choose][-1], flux[choose][-1], 
            0, -flux[choose][-1]/2, head_width=freq[choose][-1]/6,
            length_includes_head=True, head_length=flux[choose][-1]/5,
            color='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=131d", transform=ax.transAxes)


def plot_sep_panels():
    fig,axarr = plt.subplots(4, 1, figsize=(5,8), sharex=True, sharey=True)
    islim, tel, freq, days, flux, eflux = get_data_all()
    day24(axarr[0])
    day36(axarr[1])
    day46(days, freq, flux, eflux, axarr[2])
    day71(days, freq, flux, eflux, axarr[3])
    #day94(days, freq, flux, eflux, axarr[3,0])
    #day131(days, freq, flux, eflux, axarr[3,1])

    # Formatting
    for ax in axarr:
        ax.set_ylabel("Flux Density [mJy]")
    ax = axarr[-1]
    ax.set_xlabel("Frequency [GHz]")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.01, 1.8)
    ax.set_xlim(3, 230)

    # Final formatting
    plt.tight_layout()
    plt.show()
    #plt.savefig("radio_sed_evolution.png", dpi=300)
    #plt.close()



