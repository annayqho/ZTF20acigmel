""" Peak lum phase space from Soderberg (2009) with AT2018cow
and Swift J1644 on top

Combine with the Chevalier diagram

These are two classic figures
"""

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15

squaresize = 50

def vele(ax):
    # only use this to plot the GRBs
    direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    inputf = direc + "/Soderberg2009Fig4.1.txt"

    dat = Table.read(inputf, format='ascii')
    x = dat['col1']
    y = dat['col2']

    choose = np.logical_and(x > 1, y > 3E49)
    #ax.scatter(
    #        x[choose], y[choose], marker='x', c='k', s=50, label="GRBs")

    # Swift TDE
    # just use one epoch
    # epoch 3: day 18

    # ASASSN14li
    # Using Day 143, the first day the peak is resolved
    v = 0.05
    E = 9.385748112535813e+47
    ax.scatter(
            v, E, marker='o', edgecolor='k', facecolor='k', s=50, label=None)
    ax.text(
            v/1.5, E/1.1, "ASASSN14li", fontsize=12,
            horizontalalignment='center',
            verticalalignment='top')
 
    # SN 2003L
    v = 0.12
    E = 6.969524134247872e+47
    ax.scatter(
            v, E, marker='+', c='k', s=100, label="SNe Ibc")
    ax.text(
            v*1.1, E, "2003L", fontsize=12,
            horizontalalignment='left',
            verticalalignment='center')

    # SN 2007bg
    v = 0.19
    E = 2.4548568053628734e+48
    ax.scatter(
            v, E, marker='+', c='k', s=100, label=None)
    ax.text(
            v*1.3, E, "2007bg", fontsize=12,
            verticalalignment='center',
            horizontalalignment='left')

    # SN 2003bg
    v = 0.11
    E = 8.736630286056538e+47
    ax.scatter(
            v, E, marker='+', c='k', s=100, label=None)
    ax.text(
            v, E*1.1, "2003bg", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 1998bw
    v = 1.26
    E = 4.7885240374184584e+48
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label="LLGRB-SNe")
    ax.text(
            v*1.2, E, "1998bw", fontsize=12,
            verticalalignment='center',
            horizontalalignment='left')

    # SN 2009bb
    v = 0.71
    E = 2.9571483301386266e+48
    ax.scatter(
            v, E, marker='+', c='k', s=100)
    ax.text(
            v*1.2, E, "2009bb", fontsize=12,
            verticalalignment='center',
            horizontalalignment='left')

    # SN 2006aj
    v = 2.12
    E = 7.463410905948252e+47
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            v, E*1.1, "2006aj", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2010bh
    v = 0.33
    E = 9.092980420991635e+47
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            v*1.2, E, "2010bh", fontsize=12,
            verticalalignment='center',
            horizontalalignment='left')


    # SN 1988Z
    v = 0.011
    E = 1.985096671996441e+48
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='white', label="SNe II")
    ax.text(
            v, E*1.1, "88Z", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 1979C
    v = 0.016
    E = 1.0369335161200779e+48
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='white', label=None)
    ax.text(
            v, E*1.1, "79C", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')


    # AT2018cow
    v = 0.126807067696
    E = 3.68387786116342e+48
    ax.scatter(
            v, E, 
            marker='*', s=300, facecolors='black', edgecolors='black')
    ax.text(
            v, E*1.1, "AT2018cow", 
            fontsize=14, 
            verticalalignment='bottom', horizontalalignment='center')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(
        "Blastwave Velocity $(\\Gamma \\beta)$",
        fontsize=14)
    ax.set_ylabel(
            "Energy (erg) $= U_B/\epsilon_B$, \qquad $\epsilon_e=\epsilon_B=0.33$", 
            fontsize=14)
    ax.set_ylim(5E47, 1.3E49)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(loc='upper left', fontsize=10)


    # make a twin axis
    ax2 = ax.twinx()
    ax2.set_ylabel(
            "Energy (erg) $= U_B/\epsilon_B$, $\qquad \epsilon_e=0.1;\epsilon_B=0.01$", 
            fontsize=14, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i*9
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=14)
    ax2.set_xlim(4E-3, 8)
    ax.set_xlim(4E-3, 8)

fig,ax = plt.subplots(1,1, figsize=(5,6), dpi=100)
vele(ax)

#plt.show()
plt.savefig("vel_e.eps", format='eps', bbox_inches='tight')
plt.close()
