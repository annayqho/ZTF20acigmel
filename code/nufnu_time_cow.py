""" 
Plot of luminosity over time
"""


import matplotlib
from matplotlib import rc
rc("font", family="serif")
#rc("text", usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/astronomy/papers_complete/AT2018cow/code")
sys.path.append("/Users/annaho/Dropbox/astronomy/projects_proto/IcBL/old/data/radio_compilations/Zauderer2011/")
from astropy.cosmology import Planck15
from get_radio import *
from scale_fluxes import sma_lc
from read_table import *


def at2018cow(ax, col, legend):
    """ Close to 100 GHz values """
    d = Planck15.luminosity_distance(z=0.014).cgs.value
    t = 22
    f = 94E-3*1E-23
    l = f*4*np.pi*d**2
    ax.scatter(t, l, c=col, label=legend)

    


def tde(ax, col, legend):
    """  Plot the 225/230 GHz light curve from the SMA
    
    Plot the 4.9 GHz light curve from the VLA
    """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value

    # In the Eftekhari paper, it says that although the event was first
    # trigge by Swift/BAT on 2011 March 28.55 UT, subsequent
    # analysis of the BAT data revealed discernible emission as early as
    # 2011 March 25. All times should therefore be shifted relative to Mar 25.5

    # Need to add 3.04 to the Zauderer points
    nu, dt, f, ef, islim = zauderer()
    t = (dt+3.04)/(1+z)

    nu_plt = 225E9
    choose = np.logical_and(~islim, nu == nu_plt/1E9)
    dt_all = t[choose]
    nufnu_all = f[choose]*nu_plt

    # Berger2012
    nu_plt = 230E9
    t_plt = (np.array([10.30, 11.13, 17.23, 18.25, 20.24, 21.25, 125.05]))/(1+z)
    f_plt = np.array([14.90, 11.70, 13.30, 9.90, 8.20, 8.30, 6.10])
    dt_all = np.append(dt_all, t_plt)
    nufnu_all = np.append(nufnu_all, f_plt*nu_plt)

    order = np.argsort(dt_all)
    lum = plot_line(
            ax, d, dt_all[order], nufnu_all[order], 
            'SwiftJ1644+57', 'TDE', col, legend)
    ax.text(dt_all[order][-1]*1.1, lum[-1], 'Swift J1644+57', fontsize=10,
            verticalalignment='center',
            horizontalalignment='left')

    # Low frequency
    nu_plt = 4.9E9
    choose = np.logical_and(~islim, nu == nu_plt/1E9)
    dt_all = t[choose]
    nufnu_all = nu_plt*f[choose]

    # adding the set from Berger2012
    # and making the same correction as above
    t = (np.array([3.87, 4.76, 5.00, 5.79, 6.78, 7.77, 9.79, 14.98, 22.78,
        35.86, 50.65, 67.61, 94.64, 111.62, 126.51, 143.62, 164.38, 174.47,
        197.41, 213.32])) / (1+z)
    f = np.array([0.25, 0.34, 0.34, 0.61, 0.82, 1.48, 1.47, 1.80, 2.10, 4.62,
        4.84, 5.86, 9.06, 9.10, 9.10, 11.71, 12.93, 12.83, 13.29, 12.43])
    dt_all = np.append(dt_all, t)
    nufnu_all = np.append(nufnu_all, f*nu_plt)

    # adding the set from Zauderer2013
    # they also say it's relative to March 25.5...
    # so I think I need to subtract 3.04 days from here too
    t = (np.array([245.23, 302.95, 383.92, 453.66, 582.31]))/(1+z)
    f = np.array([12.17, 12.05, 12.24, 11.12, 8.90])
    dt_all = np.append(dt_all, t)
    nufnu_all = np.append(nufnu_all, f*nu_plt)

    # adding the set from Eftekhari 2018
    t = np.array([645, 651.1, 787.6, 1032, 1105, 1373, 1894])
    f = np.array([8.24, 8.63, 6.23, 4.21, 3.52, 2.34, 1.47])
    dt_all = np.append(dt_all, t)
    nufnu_all = np.append(nufnu_all, f*nu_plt)

    order = np.argsort(dt_all)
    lum = plot_line(
            ax, d, dt_all[order], nufnu_all[order], 
            'SwiftJ1644+57', 'TDE', col, legend)
    ax.text(dt_all[order][0], lum[0]/2, 'Swift J1644+57', fontsize=10,
            verticalalignment='top',
            horizontalalignment='center')


def asassn14li(ax, col, legend):
    """ Alexander et al. 2016 """
    nu = 5.0E9
    d = Planck15.luminosity_distance(z=0.0206).cgs.value
    t = np.array([80, 141.38, 207.33, 246.25, 303.01, 375.94, 389.96])
    flux = np.array([2, 1.91, 1.74, 1.56, 1.26, 0.81, 0.89])
    lum = plot_line(
            ax, d, t, nu*flux, 'ASASSN14li', 'TDE', col, legend,
            zorder=10)
    ax.text(t[-1]/2, lum[-1]/1.5, 'ASASSN14li', fontsize=10,
            verticalalignment='top',
            horizontalalignment='left')



def sn1993J(ax, col, legend):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    freq = 99.4E9
    d = 1.1E25
    nu, dt, f, ef, islim = read_1993J_high_freq()
    l = f*1E-3*1E-23*4*np.pi*d**2
    choose = np.logical_and(~islim, nu==freq)
    ax.plot(dt[choose], l[choose], c=col, label=legend, ls='--')


def sn2011dh(ax, col, legend):
    """ SN 2011dh
    Horesh et al. 2013
    Krauss et al. 2012
    M51: d = 8.03 Mpc; expl date May 31.58
    """
    d = 2.5E25

    # HIGH FREQUENCY
    # use two freq: 107E9 and 93E9
    dt, nu, f, ef, islim = read_2011dh()
    choose = np.logical_and(~islim, np.logical_or(nu==107E9, nu==93E9))
    l = f*1E-3*1E-23*4*np.pi*d**2
    ax.plot(dt[choose], l[choose], c=col, ls='--')


def sn2020oi(ax, col, legend):
    """ Maeda et al. 2021 """
    d = 15.5 * 3.086E24
    dt = np.array([5.4, 8.4, 18.3, 51.3])
    fnu = np.array([1.3, 1.22, 0.196, 0.115])
    l = fnu*1E-3*1E-23*4*np.pi*d**2
    ax.plot(dt, l, c=col, ls='--')


def grb030329(ax, col, legend):
    """ Sheth et al. 2003 
    Berger 2003
    Van der Horst et al. 2008
    
    Explosion day was obviously 03/29
    """
    z = 0.1686
    d = Planck15.luminosity_distance(z=z).cgs.value

    # HIGH FREQUENCY

    # Sheth
    dat = np.loadtxt("030329_dat.txt", delimiter=',', dtype=str)
    freq = 100E9
    t = dat[:,0].astype(float)
    lum = dat[:,1].astype(float) * 1E-3 * 1E-23 * 4 * np.pi * d**2
    ax.plot(t, lum, c=col)


def grb181201A(ax, col, legend):
    """ Laskar et al. 2018
    """
    z = 0.450
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 97
    t = np.array([8.847e-01, 1.917e+00, 3.877e+00, 8.661e+00, 2.981e+01])
    flux = np.array([3.413, 1.987, 1.199, 0.624, 0.259])
    lum = flux * 1E-3 * 1E-23 * 4* np.pi*d**2
    ax.plot(t, lum, c=col)


def grb161219B(ax, col, legend):
    """ Laskar et al. """
    z = 0.1475
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 100
    t = np.array([1.30, 3.30, 8.31, 24.45, 78.18])
    flux = np.array([1244, 897, 500, 285, 51])
    lum = flux * 1E-6 * 1E-23 * 4* np.pi*d**2
    ax.plot(t, lum, c=col)



def grb130427A(ax, col, legend):
    """ Perley et al
    They have data from CARMA/PdBI at 90 GHz (3mm)
    and also CARMA at 85.00 GHz, which is close enough
    But by the time they caught it, it was fading
    """
    z = 0.340
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 93E9
    obs_t_1 = np.array([0.77, 1, 1.91, 2.8]) 
    obs_flux_1 = np.array([3416, 2470, 1189, 807]) * 1E-6
    obs_lum_1 = obs_flux_1 * 1E-23 * 4 * np.pi * d**2

    freq = 85E9 
    obs_t_2 = np.array([0.81, 3.58, 6.41, 10.36, 23.52])
    obs_flux_2 = np.array([3000, 903, 588, 368, 197]) * 1E-6
    obs_lum_2 = obs_flux_2 * 1E-23 * 4 * np.pi * d**2

    t = np.hstack((obs_t_1, obs_t_2))
    lum = np.hstack((obs_lum_1, obs_lum_2))
    order =np.argsort(t)

    ax.plot(obs_t_2, obs_lum_2, c=col, label=legend)


def sn1998bw(ax, col, legend):
    """ SN 1998bw
    
    This is a bit complicated because there are two peaks in the light curve,
    but I am using a frequency that has a main peak rather than a frequecy
    with two clear distinct peaks...
    """
    d = 1.17E26 # cm
    nu = 150E9
    t = np.array([12.4])
    f = np.array([39])
    lum = f*1E-3*1E-23*4*np.pi*d**2
    ax.scatter(t, lum, c=col, label=legend)



def sn2017iuk(ax, col, legend):
    """ SN 2017iuk
    """
    d = Planck15.luminosity_distance(z=0.0368).cgs.value
    nu = 92E9 # Band 3
    t = np.array([6.10])
    f = np.array([28])
    l = f*1E-3*1E-23*4*np.pi*d**2
    ax.scatter(t,l,c=col)


def at2020xnd(ax, col, legend):
    dcm = Planck15.luminosity_distance(z=0.2442).cgs.value
    dt = [17,24,31,39,46,67]
    fnu = np.array([305,648,1030,868,558,330])
    efnu = np.array([57,44,44,46,38,48])
    nu = np.array([94,94,94,94,94,79])
    lum = fnu * 1E-6 * 1E-23 * 4 * np.pi * dcm**2 
    ax.plot(dt, lum, c=col, ls=':', lw=4, label=legend)


if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(5,7), dpi=100)
    props = dict(boxstyle='round', facecolor='white')

    sn_col = '#7a5195'
    llgrb_col = '#ef5675'
    cow_col = '#ffa600'
    lgrb_col = '#003f5c'

    # First category: long-duration GRBs
    grb130427A(ax, lgrb_col, legend='LGRB')
    grb030329(ax, lgrb_col, None)
    grb181201A(ax, lgrb_col, None)
    grb161219B(ax, lgrb_col, None)

    # Second category: low-luminosity GRBs
    sn1998bw(ax, llgrb_col, legend='LLGRB')
    sn2017iuk(ax, llgrb_col, None)

    # Third category: Cow-like
    #at2018cow(ax, cow_col, 'Cow-like')
    at2020xnd(ax, cow_col, 'AT2020xnd')

    # Final category: 'ordinary' CC SNe
    sn1993J(ax, sn_col, 'CC SN')
    sn2011dh(ax, sn_col, None)
    sn2020oi(ax, sn_col, None)

    ax.set_ylabel(
            r"Luminosity $L_{\nu}$ [erg$\,$s$^{-1}$Hz$^{-1}$]", 
            fontsize=16)
    ax.set_title(r"$\nu \approx 90\,\mathrm{GHz}$", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlim(0.7, 300) 
    ax.set_ylim(1E25, 2E32)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Time [days; observer frame]", fontsize=16)
    ax.legend(fontsize=12, loc='upper right')

    plt.tight_layout()
    #plt.show()
    plt.savefig(
            "lum_evolution.png", dpi=300, 
            bbox_inches='tight', pad_inches=0.1)
    plt.close()
