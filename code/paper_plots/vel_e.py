""" Peak lum phase space from Soderberg (2009) with AT2018cow
and Swift J1644 on top

Combine with the Chevalier diagram

These are two classic figures
"""

import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Dropbox/astronomy/papers_active/ZTF20acigmel/code")
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
from format import *

d = get_format()

squaresize = 50

def vele(ax):
    # Dillon's MAXI source...Dong et al. 2021
    v = 0.025
    E = 2E49
    ax.scatter(
            v, E, marker='X', edgecolor='k', facecolor='k', s=50, label=None)
    ax.text(
            v/1.5, E/1.1, "VT1210+4956", fontsize=d['font_small'],
            horizontalalignment='center',
            verticalalignment='top')

    # FIRST transient
    # Mooley+18, Law+2021
    v = 0.02
    E = 3E49
    ax.scatter(
            v, E, marker='X', edgecolor='k', facecolor='k', s=50, label=None)
    ax.text(
            v/1.5, E/1.1, "FIRST J1419", fontsize=d['font_small'],
            horizontalalignment='center',
            verticalalignment='top')

    # ASASSN14li
    # Using Day 143, the first day the peak is resolved
    v = 0.05
    E = 9.385748112535813e+47
    ax.scatter(
            v, E, marker='o', edgecolor='k', facecolor='k', s=50, label=None)
    ax.text(
            v/1.5, E/1.1, "ASASSN14li", fontsize=d['font_small'],
            horizontalalignment='center',
            verticalalignment='top')

    # PTF11qcj
    v = 0.05
    E = 3E48
    ax.scatter(
            v, E, marker='+', edgecolor='k', facecolor='k', s=100, label=None)
    ax.text(
            v/1.5, E/1.1, "PTF11qcj", fontsize=d['font_small'],
            horizontalalignment='center',
            verticalalignment='top')
 
    # SN 2003L
    v = 0.12
    E = 6.969524134247872e+47
    ax.scatter(
            v, E, marker='+', c='k', s=100, label="SNe Ibc")
    ax.text(
            v*1.1, E, "2003L", fontsize=d['font_small'],
            horizontalalignment='left',
            verticalalignment='center')

    # SN 2007bg
    v = 0.19
    E = 2.4548568053628734e+48
    ax.scatter(
            v, E, marker='+', c='k', s=100, label=None)
    ax.text(
            v*1.3, E, "2007bg", fontsize=d['font_small'],
            verticalalignment='center',
            horizontalalignment='left')

    # SN 2003bg
    v = 0.11
    E = 8.736630286056538e+47
    ax.scatter(
            v, E, marker='+', c='k', s=100, label=None)
    ax.text(
            v, E*1.1, "2003bg", fontsize=d['font_small'],
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 1998bw
    v = 1.26
    E = 4.7885240374184584e+48
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label="LLGRB-SNe")
    ax.text(
            v*1.2, E, "1998bw", fontsize=d['font_small'],
            verticalalignment='center',
            horizontalalignment='left')

    # SN 2009bb
    v = 0.71
    E = 2.9571483301386266e+48
    ax.scatter(
            v, E, marker='+', c='k', s=100)
    ax.text(
            v*1.2, E, "2009bb", fontsize=d['font_small'],
            verticalalignment='center',
            horizontalalignment='left')

    # SN 2006aj
    v = 2.12
    E = 7.463410905948252e+47
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            v, E*1.1, "2006aj", fontsize=d['font_small'],
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2010bh
    v = 0.33
    E = 9.092980420991635e+47
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            v*1.2, E, "2010bh", fontsize=d['font_small'],
            verticalalignment='center',
            horizontalalignment='left')


    # SN 1988Z
    v = 0.011
    E = 1.985096671996441e+48
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='white', label="SNe II")
    ax.text(
            v, E*1.1, "88Z", fontsize=d['font_small'],
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 1979C
    v = 0.016
    E = 1.0369335161200779e+48
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='white', label=None)
    ax.text(
            v, E*1.1, "79C", fontsize=d['font_small'],
            verticalalignment='bottom',
            horizontalalignment='center')


    # AT2018cow
    v = 0.126807067696
    E = 3.68387786116342e+48
    col = d['colors']['4'][1]
    ax.scatter(
            v, E, 
            marker='o', s=100, c=col)
    ax.text(
            v, E*1.1, "AT2018cow", color=col,
            fontsize=d['font_small'], 
            verticalalignment='bottom', horizontalalignment='center')

    v = 0.22849036502452907
    E = 6.2E48
    col = d['colors']['4'][2]
    ax.scatter(
            v, E, marker='*', s=300, c=col)
    ax.text(
            v, E*1.1, "AT2020xnd", 
            fontsize=d['font_small'], color=col,
            verticalalignment='bottom', horizontalalignment='center')

    # Koala
    v = 0.4
    E = 3E49
    col = d['colors']['4'][0]
    ax.scatter(
            v, E, marker='D', s=100, c=col)
    ax.text(
            v, E*1.1, "Koala", 
            fontsize=d['font_small'], color=col,
            verticalalignment='bottom', horizontalalignment='center')

    v = 0.55
    E = 5.6E49
    col = d['colors']['4'][3]
    ax.scatter(
            v, E, marker='h', s=100, c=col)
    ax.text(
            v, E*1.1, "CSS161010", 
            fontsize=d['font_small'], color=col,
            verticalalignment='bottom', horizontalalignment='center')


    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(
        "Blastwave Velocity $(\\Gamma \\beta)$",
        fontsize=d['font_med'])
    ax.set_ylabel(
            "Energy (erg) $= U_B/\epsilon_B$, $\qquad \epsilon_e=\epsilon_B=0.33$", 
            fontsize=d['font_med'])
    ax.set_ylim(5E47, 9E49)
    ax.tick_params(axis='both', labelsize=d['font_med'])
    ax.legend(loc='upper left', fontsize=d['font_small'])


    # make a twin axis
    ax2 = ax.twinx()
    ax2.set_ylabel(
            r"Energy (erg) $= U_B/\epsilon_B$, $\qquad \epsilon_e=0.1;\epsilon_B=0.01$", 
            fontsize=d['font_med'], rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i*9
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=d['font_med'])
    ax2.set_xlim(4E-3, 8)
    ax.set_xlim(4E-3, 8)

fig,ax = plt.subplots(1,1, figsize=(5,5), dpi=100)
vele(ax)

plt.tight_layout()
#plt.show()
plt.savefig("vel_e.png", dpi=200)
plt.close()
