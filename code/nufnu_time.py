""" 
Plot of luminosity over time
"""


import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/astronomy/papers_complete/AT2018cow/code")
sys.path.append("/Users/annaho/Dropbox/astronomy/papers_complete/AT2018cow/data")
sys.path.append("/Users/annaho/Dropbox/astronomy/projects_proto/IcBL/data/radio_compilations/Zauderer2011")
from astropy.cosmology import Planck15
from get_radio import *
from scale_fluxes import sma_lc
from read_table import *


def plot_limits(ax, x, y, ratiox, ratioy, col):
    """ Plot two arrows from the point """
    ax.annotate('', xy=(x*ratiox, y), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(x, y*ratioy), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))


def plot_line(ax, d, t, nufnu, name, label, col, legend=False, zorder=1):
    """ Plot a line
    If nu > 90 GHz, make it on the left axis
    If nu < 10 GHz, make it on the right axis
    
    Parameters
    ----------
    nufnu: Hz * mJy 
    name: name of the source
    label: label to use as the legend (which also determines the col)
    """
    lum = nufnu * 1e-23 * 1e-3 * 4 * np.pi * d**2
    fs = 11
    nsize = 10 # normal size for points
    if name=='AT2018cow':
        marker='*'
        fcol = col
        s=70
    else:
        if label=='SN':
            marker='o'
            s=nsize
            fcol = col # fill color
            label = 'Supernova'
        elif label=='GRB':
            marker='o'
            fcol = 'white' # unfilled
            s=nsize
            label = 'Gamma-ray burst'
        elif label=='Rel. SN':
            marker='s'
            fcol = col 
            s=nsize
            label = 'Relativistic supernova'
        elif label=='TDE':
            marker='s'
            fcol = 'white' #unfilled
            s=nsize
            label='Tidal disruption event'
    ax.scatter(
            t, lum, facecolor=fcol, edgecolor=col, 
            marker=marker, s=s, zorder=zorder)
    if legend:
        ax.plot(t, lum, c=col, ls='-', label=label, zorder=zorder)
    else:
        ax.plot(t, lum, c=col, ls='-', label=None, zorder=zorder)
    return lum


def plot_points(ax, d, nu, t, f, marker, name=None):
    """ Plot set of two points """
    lums = []
    for ii,nuval in enumerate(nu):
        if nuval > 90E9:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker, name=name)
            lums.append(lum)
        else:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker)
            lums.append(lum)
    ax.plot(
        t, lums, ls='--', c='k')
    return lums


def koala(ax, col, legend):
    # This is the X-band light curve
    dt = np.array([81,310,352,396])
    f = np.array([8.3E39,1.5E39,1.1E39,8.16E38])-2.65E38
    ax.errorbar(dt/(1.2714), f, 0.0006*f, c='black', fmt='*', ms=10)
    ax.plot(dt/(1.2714), f, c='black', lw=2)
    ax.text(
            dt[-1]*1.5, f[-1]/1.5, "ZTF18abvkwla",
            fontsize=10, horizontalalignment='center',
            verticalalignment='top')


def at2018cow(ax, col, legend):
    """ 231.5 GHz light curve and 9 GHz light curve """
    d = Planck15.luminosity_distance(z=0.014).cgs.value

    # high frequency
    a, b, c = sma_lc()
    dt, f, ef = b
    ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
    nu = 231.5E9
    lum = plot_line(
            ax[0], d, dt, nu*f, 
            'AT2018cow', None, col, legend, zorder=10)
    ax[0].text(dt[0]/1.2, lum[0], 'AT2018cow', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    # low frequency
    nu = 9E9
    dat = Table.read(
        "/Users/annaho/Dropbox/astronomy/papers_complete/AT2018cow/data/radio_lc.dat", delimiter="&",
        format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = np.logical_or(tel == 'SMA', tel == 'ATCA')

    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    eflux_sys = np.array([0.1*f for f in flux])
    eflux_form = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])
    eflux = np.sqrt(eflux_sys**2 + eflux_form**2)
    choose = freq == 9


    # add the Margutti point and the Bietenholz point
    margutti_x = np.array([84,287])
    margutti_y = np.array([6E28, 3.2E26])/(4*np.pi*d**2)/1E-23/1E-3
    x = np.hstack((days[choose], margutti_x))
    y = np.hstack((flux[choose], margutti_y)) * nu
    lum = plot_line(
            axarr[1], d, x, y,
            'AT2018cow', None, col, legend, zorder=10)


def maxi(ax):
    """ data from Dillon """
    d = 4.62E26
    t = 4*365
    nu = 6E9
    f = 6.16
    lum = plot_point(ax, d, nu, t, f, '*')
    plt.text(t/1.2, lum, 'MAXI 140814A', fontsize=11,
            horizontalalignment='right', verticalalignment='center')


def tde(ax, col, legend):
    """  Plot the 225/230 GHz light curve from the SMA
    
    Plot the 4.9 GHz light curve from the VLA
    """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value

    # In the Eftekhari paper, it says that although the event was first
    # triggered by Swift/BAT on 2011 March 28.55 UT, subsequent
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
            ax[0], d, dt_all[order], nufnu_all[order], 
            'SwiftJ1644+57', 'TDE', col, legend)
    ax[0].text(dt_all[order][-1]*1.1, lum[-1], 'Swift J1644+57', fontsize=11,
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
            ax[1], d, dt_all[order], nufnu_all[order], 
            'SwiftJ1644+57', 'TDE', col, legend)
    ax[1].text(dt_all[order][0], lum[0]/2, 'Swift J1644+57', fontsize=11,
            verticalalignment='top',
            horizontalalignment='center')


def asassn14li(ax, col, legend):
    """ Alexander et al. 2016 """
    nu = 5.0E9
    d = Planck15.luminosity_distance(z=0.0206).cgs.value
    t = np.array([80, 141.38, 207.33, 246.25, 303.01, 375.94, 389.96])
    flux = np.array([2, 1.91, 1.74, 1.56, 1.26, 0.81, 0.89])
    lum = plot_line(
            ax[1], d, t, nu*flux, 'ASASSN14li', 'TDE', col, legend,
            zorder=10)
    ax[1].text(t[-1]/2, lum[-1]/1.5, 'ASASSN14li', fontsize=11,
            verticalalignment='top',
            horizontalalignment='left')


def sn2003L(ax, col, legend):
    """ Soderberg et al
    Values at 8.5 GHz """
    d = 2.8432575937224894e+26
    nu_plt = 8.5E9
    nu, dt, f, ef = read_2003L()
    choose = nu == nu_plt
    lum = plot_line(
            ax[1], d, dt[choose], 1E-3*f[choose]*nu_plt, 
            'SN2003L', 'SN', col, legend)
    #ax[1].text(dt[choose][-1]/1.05, lum[-1], 'SN2003L', fontsize=11,
    #        verticalalignment='center',
    #        horizontalalignment='left')
    

def sn1979c(ax, col, legend):
    """ Weiler 1986 and Weiler 1991
    This is a IIL 
    
    Too old to have the raw table on arXiv
    """
    d = 5.341805643483106e+25
    nu = 1.4E9 # 20cm
    t = np.array(
            [437,594,631,663,679,684,727,747,786,822,839,876,882,
             914,937,973,995,
             1026,1071,1091,1127,1156,1168,1212,1243,1277,1314,1358,1390,
             1415,1435,1466,1513,1565,1600,1634,1659,1698,1714,1750,1771,
             1931,2027])
    flux = np.array(
            [0.2,2.1,2.5,2.7,2.8,2.8,4.4,4.8,6.0,7.1,7.1,7.6,8.6,
             9.8,6.5,8.6,9.5,
             10.2,10.8,10.3,10.4,12.2,10.1,10.2,11.5,11.2,13.0,11.3,10.2,
             9.6,11.2,13.2,11.1,9.1,8.5,9.1,8.8,10.1,9.7,9.1,8.9,
             7.0,7.7])
    lum = plot_line(ax[1], d, t, nu*flux, 'SN1979c', 'SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], 'SN1979C', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')
    

def sn1993J(ax, col, legend):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    d = 1.1e25
    freq = 5E9
    nu, dt, f, ef, islim = read_1993J_low_freq()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[1], d, dt[choose], freq*f[choose], 
            'SN1993J', 'SN', col, legend)
    ax[1].text(dt[choose][0]/1.05, lum[0], 'SN1993J', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    freq = 99.4E9
    nu, dt, f, ef, islim = read_1993J_high_freq()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[0], d, dt[choose], freq*f[choose], 
            'SN1993J', 'SN', col, legend)
    ax[0].text(dt[choose][-1]*1.1, lum[-1], 'SN1993J', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')


def sn2011dh(ax, col, legend):
    """ SN 2011dh
    Horesh et al. 2013
    Krauss et al. 2012
    M51: d = 8.03 Mpc; expl date May 31.58
    """
    d = 2.5E25

    dt, nu, f, ef, islim = read_2011dh()

    # HIGH FREQUENCY
    # use two freq: 107E9 and 93E9
    choose = np.logical_and(~islim, np.logical_or(nu==107E9, nu==93E9))
    lum = plot_line(
            ax[0], d, dt[choose], nu[choose]*f[choose], 
            'SN2011dh', 'SN', col, legend)
    ax[0].text(dt[choose][0]/1.05, lum[0], 'SN2011dh', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    # LOW FREQUENCY
    # from Horesh 2013 and Krauss 2012
    dt_all = []
    f_all = []
    nu_all = []

    freq = 8.5E9
    choose = np.logical_and(~islim, nu==freq)
    dt_all.extend(dt[choose])
    f_all.extend(f[choose])
    nu_all.extend([freq]*sum(choose))

    freq = 6.7E9
    dt_all.extend([16.4, 20.4, 25.4, 35.3, 45.3, 58.2, 92.9])
    f_all.extend([4.09, 4.8, 5.98, 7.222, 6.987, 6.11, 3.941])
    nu_all.extend([freq]*7)

    dt_all = np.array(dt_all)
    f_all = np.array(f_all)
    nu_all = np.array(nu_all)

    lum = plot_line(
            ax[1], d, dt_all, nu_all*f_all, 
            'SN2011dh', 'SN', col, legend)
    ax[1].text(dt[choose][0]/1.05, lum[0], 'SN2011dh', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


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
    freq = 250E9
    t = np.array([1, 4, 6, 7, 8, 16, 22]) / (1+z)
    f = np.array([49.2, 45.7, 41.6, 32, 25.5, 9.3, 5.2])
    lum = plot_line(ax[0], d, t, freq*f, 'GRB030329', 'GRB', col, legend)
    ax[0].text(t[-1]*1.05, lum[-1], ' GRB 030329', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')

    # LOW FREQUENCY

    # Berger: this is the best frequency to pick from this paper
    t = np.array(
            [0.58, 1.05, 2.65, 3.57, 4.76, 6.89, 7.68, 9.49, 11.90, 
                12.69, 14.87, 16.66, 18.72, 20.58, 25.70, 28.44, 31.51, 
                33.58, 36.52, 42.55, 44.55, 59.55, 66.53]) / (1+z)
    f = np.array(
            [3.50, 1.98, 8.50, 6.11, 9.68, 15.56, 12.55, 13.58, 17.70, 
                17.28, 19.15, 17.77, 15.92, 16.08, 15.34, 12.67, 13.55, 
                13.10, 10.64, 8.04, 8.68, 4.48, 4.92])
    nu = np.array([8.5E9]*len(f))

    # Van der Horst: best frequency is 2.3 GHz
    t = np.append(t, np.array([268.577, 306.753, 365.524, 420.168, 462.078, 
        583.683, 743.892, 984.163]) / (1+z))
    f = np.append(
            f, np.array([1613, 1389, 871, 933, 707, 543, 504, 318]) * 1E-3)
    nu = np.append(nu, np.array([2.3E9]*8))
    lum = plot_line(ax[1], d, t, nu*f, 'GRB030329', 'GRB', col, legend)
    ax[1].text(t[-8]*1.05, lum[-8]/1.1, 'GRB030329', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='left')
    


def grb130427A(ax, col, legend):
    """ Perley et al
    They have data from CARMA/PdBI at 90 GHz (3mm)
    But by the time they caught it, it was fading
    """
    z = 0.340
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 93E9
    t = np.array([0.77, 1, 1.91, 2.8]) / (1+z)
    flux = np.array([3416, 2470, 1189, 807]) * 1E-3
    lum = plot_line(ax[0], d, t, freq*flux, 'GRB130427A', 'GRB', col, legend)
    ax[0].text(t[-1]*1.05, lum[-1], 'GRB130427A', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')

    freq = 5.10E9
    t = np.array([0.677, 2.04, 4.75, 9.71, 17.95, 63.78, 128.34]) / (1+z)
    f = np.array([1290, 1760, 648, 454, 263, 151, 86]) * 1E-3
    lum = plot_line(ax[1], d, t, freq*f, 'GRB130427A', 'GRB', col, legend)
    ax[1].text(t[-3]*1.05, lum[-3], 'GRB130427A', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='left')


def sn2007bg(ax, col, legend):
    """ Salas et al. 2013
    Peak is resolved for 4.86, 8.46 GHz
    """
    nu = 8.46E9
    d = Planck15.luminosity_distance(z=0.0346).cgs.value
    t = np.array(
            [13.8, 19.2, 26.1, 30.9, 41.3, 55.9, 66.8, 81.8, 98.8, 124, 
                144, 159.8, 189.9, 214.9, 250.9, 286.8, 314.8, 368.8, 
                386.8, 419.9, 566.9, 623.8, 720.8, 775.8, 863.8])
    f = np.array(
            [480, 753, 804, 728, 1257, 1490, 1390, 1325, 1131, 957, 
                621, 316, 379, 404, 783, 1669, 2097, 2200, 
                2852, 3344, 3897, 3891, 3842, 3641, 3408]) * 1E-3
    lum = plot_line(ax[1], d, t, nu*f, 'SN2007bg', 'SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], 'SN2007bg', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right')


def sn2003bg(ax, col, legend):
    """ Soderberg et al.
    Peak is resolved for 22.5, 15, 8.46, 4.86, 1.43
    Again, there are two peaks...
    Let's choose the first peak, 8.46
    """
    nu = 8.46E9
    d = 6.056450393620008e+25

    t = np.array(
                [10, 12, 23, 35, 48, 58, 63, 73, 85, 91, 115, 129,
                132, 142, 157, 161, 181, 201, 214, 227, 242, 255,
                266, 285, 300, 326, 337, 351, 368, 405, 410, 424,
                434, 435, 493, 533, 632, 702, 756, 820, 902, 978])
    f = np.array(
                [2.51, 3.86, 12.19, 24.72, 40.34, 51.72, 49.64, 46.20,
                38.638, 33.85, 45.74, 53.94, 54.27, 54.83, 48.43,
                47.43, 35.76, 31.35, 28.67, 27.38, 24.57, 22.30,
                21.67, 21.31, 20.88, 20.33, 19.85, 18.84, 17.14,
                14.61, 14.49, 14.16, 13.25, 13.08, 10.04, 8.92,
                6.23, 6.18, 4.62, 3.93, 4.69, 4.48])
    lum = plot_line(ax[1], d, t, nu*f, 'SN2003bg', 'SN', col, legend)
    #ax[1].text(t[0]/1.05, lum[0], 'SN2003bg', fontsize=11,
    #        verticalalignment='center',
    #        horizontalalignment='right')


def sn2009bb(ax, col, legend):
    """ expl date Mar 19 """
    nu = 8.46E9
    d = 1.237517263280789e+26
    t_apr = 11 + np.array([5.2, 8.2, 13.2, 15.1, 23.2, 29.1])
    t_may = 11 + 30 + np.array([3.1, 10.1, 13, 20.1, 27])
    t_jun = 11 + 30 + 31 + np.array([6, 17, 26])
    t_jul = 11 + 30 + 31 + 30 + np.array([18.9])
    t_aug = 11 + 30 + 31 + 30 + 31 + np.array([11.8])
    t = np.hstack((t_apr, t_may, t_jun, t_jul, t_aug))
    flux = np.array([24.681, 17.568, 16.349, 13.812, 8.881,
        7.714, 8.482, 6.824, 6.327, 3.294, 4.204, 3.203, 2.392,
        1.903, 1.032, 1.084])
    lum = plot_line(ax[1], d, t, nu*flux, 'SN2009bb', 'Rel. SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], '2009bb', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


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
    lum = plot_line(ax[0], d, t, nu*f, 'SN1998bw', 'Rel. SN', col, legend)
    nu = 2.3E9
    t = np.array([11.7, 14.6, 15.7, 16.5, 17.8, 19.7, 21.6, 23.6, 25.9, 26.8, 28.8, 30.0, 32.9, 34.7, 36.8, 38.8, 40.0, 45.7, 51.7, 57.7, 64.7, 67.7, 80.5])
    f = np.array([19.7, 22.3, 23.5, 23.9, 25.1, 25.3, 20.9, 22.9, 28.0, 28.7, 31.1, 31.3, 27.3, 33.5, 31.8, 31, 31.3, 26.8, 23.1, 18.5, 15.6, 15.6, 9.6])
    ax[0].text(t[0]*1.2, lum[0], '98bw', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')
    lum = plot_line(ax[1], d, t, nu*f, 'SN1998bw', 'Rel. SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], '1998bw', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def css(ax, col, legend):
    """ 6 GHz light curve """
    d = Planck15.luminosity_distance(z=0.034).cgs.value

    # low frequency
    nu = 6E9

    # add the points from Deanne's paper
    x = np.array([69, 99, 162, 357])
    y = np.array([4.5, 6.1, 2.3, 0.07])*nu*1E-3*1E-23*4*np.pi*d**2

    ax.errorbar(x/(1.034), y, 0.0006*y, c='black', fmt='*', ms=10)
    ax.plot(x/(1.034), y, c='black', lw=2)
    ax.text(
            x[-1]*1.9, y[-1]/1.3, "CSS161010",
            fontsize=11, horizontalalignment='center',
            verticalalignment='top')


def othersn(ax):
    """ 
    SN 2013ak
    most luminous events are typically SNe IIn (Chevalier 2006)
    SN Ib/c: Berger 2003, Chevalier & Fransson 2006
    2008D
    """




if __name__=="__main__":
    fig, axarr = plt.subplots(
            1, 2, figsize=(10,7), sharex=True, sharey=True, dpi=100)
    props = dict(boxstyle='round', facecolor='white')

    # viridis color palette
    # 440154 dark purple
    # 21918c dark green/turquoise
    # 3b528b dark blue
    # 5ec962 light green
    # fde725 yellow

    # inferno color palette
    # 000004 black
    # 420a68 dark purple
    # 932667 kind of a purple/red color
    # dd513a dark orange
    # fca50a light orange
    # fcffa4 really light yellow
    # 57106e purple
    # bc3754 pink/red
    # f98e09 orange

    # OK so I want to use four colors. I think it's nice to have two light ones
    # and two dark ones.

    #maxi(ax)
    tde(axarr, '#57106e', legend=True)
    asassn14li(axarr, '#57106e', None)

    sn2003L(axarr, 'lightblue', legend=True)
    sn1979c(axarr, 'lightblue', None)
    sn1993J(axarr, 'lightblue', None)
    sn2011dh(axarr, 'lightblue', None)
    sn2007bg(axarr, 'lightblue', None)
    sn2003bg(axarr, 'lightblue', None)

    grb030329(axarr, '#f98e09', legend=True)
    grb130427A(axarr, '#f98e09', None)

    #sn2009bb(axarr, '#bc3754', legend=True)
    sn1998bw(axarr, '#bc3754', None)

    at2018cow(axarr, 'k', None)

    # new source
    # first ZTF det was around Oct 10
    # first NOEMA obs was Oct 27
    # second NOEMA obs was Nov 3
    # third NOEMA obs was Nov 10
    dcm = Planck15.luminosity_distance(z=0.2442).cgs.value
    dt = [17,24,31,39,46]
    fnu = np.array([305,648,1030,868,558])
    efnu = np.array([57,44,44,46,38])
    lum = fnu * 1E-6 * 1E-23 * 4 * np.pi * dcm**2 * 94.245E9 
    elum = efnu * 1E-6 * 1E-23 * 4 * np.pi * dcm**2 * 94.245E9
    axarr[0].errorbar(
            dt, lum, yerr=elum, fmt='D-', mec='k', c='Crimson',
            mfc='Crimson', zorder=20, ms=6)
    axarr[0].text(
            dt[-1]*1.2, lum[-1], 'AT2020xnd', color='Crimson', fontsize=16)

    # the first VLA observation was Oct 23
    dt = [13, 25, 36, 51, 71, 94]
    fnu = np.array([21, 57, 80, 154, 158, 174])
    efnu = np.array([5, 5, 5, 5, 9, 10])
    lum = fnu * 1E-6 * 1E-23 * 4 * np.pi * dcm**2 * 8E9
    elum = efnu * 1E-6 * 1E-23 * 4 * np.pi * dcm**2 * 8E9
    axarr[1].errorbar(
            dt, lum, yerr=elum, fmt='D-', mec='k', c='Crimson',
            mfc='Crimson', zorder=20, ms=6)
    axarr[1].text(
            dt[0]/1.2, lum[0], 'AT2020xnd', color='Crimson', 
            fontsize=16, ha='right')

    # Koala
    koala(axarr[1], 'k', None)

    # CSS transient
    css(axarr[1], 'k', None)

    axarr[0].set_ylabel(
            r"Luminosity $\nu L_{\nu}$ [erg s$^{-1}$]", 
            fontsize=16)
    axarr[0].set_title(r"$\nu > 90\,\mathrm{GHz}$", fontsize=16)
    axarr[1].set_title(r"$\nu < 10\,\mathrm{GHz}$", fontsize=16)
    for ax in axarr:
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xlim(0.3, 2000) 
        ax.set_ylim(1E34, 2E43)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r"Time [days; rest frame]", fontsize=16)
    axarr[1].legend(fontsize=12, loc='upper right')

    plt.subplots_adjust(wspace=0.05)
    #plt.show()
    plt.savefig(
            "lum_evolution.png", dpi=200,
            bbox_inches='tight', pad_inches=0.1)
