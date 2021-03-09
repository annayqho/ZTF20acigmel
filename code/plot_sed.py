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
    Plot VLA and NOEMA
    """

    choose = np.logical_and(days>23, days<26)
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='o', c='k')

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
    ax.text(0.05,0.85,"$\Delta t$=24-25d", transform=ax.transAxes)


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
    choose = np.logical_and(days>35, days<40)
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='o', c='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=36-39d", transform=ax.transAxes)


def day46(days, freq, flux, eflux, ax):
    """ 
    Interpolate VLA/ATCA, fix NOEMA
    """
    choose = np.logical_and(days>45, days<52)
    ax.errorbar(freq[choose], flux[choose], eflux[choose], fmt='o', c='k')

    # SMA upper limit
    choose = days==47
    ax.arrow(freq[choose], flux[choose], 0, -flux[choose]/3, color='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=46-51d", transform=ax.transAxes)


def day71(days, freq, flux, eflux, ax):
    """ 
    Just plot VLA
    """
    day=71

    # Plot the VLA LC
    x = 8
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 12
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 17.7
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 26.5
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 36.2
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=71d", transform=ax.transAxes)


def day94(days, freq, flux, eflux, ax):
    """ 
    Just plot VLA
    """
    day=94

    # Plot the VLA LC
    x = 8
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 12
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 17.7
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 26.5
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=94d", transform=ax.transAxes)


def day131(days, freq, flux, eflux, ax):
    """ 
    Just plot VLA
    """
    day=131

    # Plot the VLA LC
    x = 4
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 8
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 12
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    x = 17.7
    choose = np.logical_and(days==day,freq==x)
    ax.errorbar(x, flux[choose], eflux[choose], marker='o', c='k')

    ax.scatter(26.5, 0.048, marker='_', c='k')
    ax.arrow(26.5, 0.048, 0, -0.015, head_length=0.005, length_includes_head=True, head_width=3, color='k')

    # Label
    ax.text(0.05,0.85,"$\Delta t$=94d", transform=ax.transAxes)


if __name__=="__main__":
    fig,axarr = plt.subplots(4, 2, figsize=(7,8), sharex=True, sharey=True)
    islim, tel, freq, days, flux, eflux = get_data_all()
    day17(days, freq, flux, eflux, axarr[0,0])
    day24(axarr[0,1])
    day31(axarr[1,0])
    day36(axarr[1,1])
    day46(days, freq, flux, eflux, axarr[2,0])
    #day71(days, freq, flux, eflux, axarr[2,0])
    #day94(days, freq, flux, eflux, axarr[2,1])
    #day131(days, freq, flux, eflux, axarr[3,0])

    # Formatting
    ax = axarr[0,0]
    for ax in axarr[:,0]:
        ax.set_ylabel("Flux Density [mJy]")
    for ax in axarr[1]:
        ax.set_xlabel("Frequency [GHz]")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.01, 1.8)
    ax.set_xlim(3, 230)

    # Final formatting
    plt.tight_layout()
    plt.show()
    #plt.savefig("radio_sed_evolution.png", dpi=300)
