""" Compare the data to other SNe: SN2003L SN2007bg SN2003bg """

import numpy as np
import matplotlib.pyplot as plt
from format import get_format

# Common dL = 1E28

# SN2003L: 92 Mpc
scale_sn2003L = 92*3.086E24 / 1E28
scale_sn2003bg = 19.6*3.086E24 / 1E28
scale_at2020xnd = 1267*3.086E24 / 1E28


def compare_lc():
    """ Compare the light curves at four well-sampled frequencies
    45 GHz / 33 GHz / 15 GHz / 10 GHz
    """
    # 4-panel figure
    fig, axarr = plt.subplots(2,2,figsize=(6,6))
    d = get_format()
    cols = d['colors']['4']

    # 45 GHz
    ax = axarr[0,0]
    ax.set_title("45 GHz")

    # AT2020xnd
    dt = np.array([36.9, 51.9, 71.9])
    f = np.array([0.675, 0.668, 0.209])
    ef = np.array([0.171, 0.168, 0.063])
    ax.errorbar(dt, f*scale_at2020xnd, ef*scale_at2020xnd, fmt='o', c=cols[2],
            label='20xnd')

    # SN2003bg
    dt = np.array([35, 48, 58, 63, 73])
    f = np.array([64, 31.75, 20.56, 15.62, 10.18])
    ef = np.array([1.95, 1.43, 5.87, 2.53, 2.55])
    #ax.errorbar(dt, f*scale_sn2003bg, ef*scale_sn2003bg, fmt='o', c=cols[0])
    ax.plot(dt, f*scale_sn2003bg, c=cols[0], label='03bg')

    # 22 GHz
    ax = axarr[0,1]
    ax.set_title("22 GHz")
    
    # AT2020xnd
    dt = np.array([71.9, 94.7, 131.6])
    f = np.array([0.484, 0.301, 0.087])
    ef = np.array([0.01, 0.065, 0.02])
    ax.errorbar(dt, f*scale_at2020xnd, ef*scale_at2020xnd, fmt='o', c=cols[2],
            zorder=20, label='20xnd')

    # SN2003L
    dt = np.array([27.3, 28.3, 29.4, 30.3, 31.2, 32.3, 41.6, 44.5, 52.3, 55.3,
        65.3, 66.4, 72.4, 78.3, 85.4, 91.1, 119, 122.1, 135.3, 147, 222.9, 506])
    f = (1/1000)*np.array([3179, 3119, 3122, 3080, 3252, 2973, 2956, 2684, 2714, 2838,
        2263, 2486, 2422, 2130, 1833, 1819, 1078, 1327, 1146, 941, 670, 280])
    ef = (1/1000)*np.array([85, 86, 75, 97, 70, 99, 96, 124, 68, 99, 73, 66, 70, 223,
        100, 141, 173, 114, 165, 161, 148, 53])
    #ax.errorbar(dt, f*scale_sn2003L, ef*scale_sn2003L, fmt='o', c=cols[1])
    ax.plot(dt, f*scale_sn2003L, c=cols[1], label='03L', ls='--')
    
    # SN2003bg
    dt = np.array([12, 23, 35, 48, 58, 63, 73, 85, 91, 115,
        129, 181, 201, 214, 227, 255])
    f = np.array([30.77, 106.30, 85.39, 74.58, 39.91, 32.48, 30.53, 9.63, 10.70,
        21.07, 30.19, 11, 6.59, 10.13, 7.84, 9.38])
    ef = np.array([0.71, 2.16, 1.72, 1.53, 0.89, 0.71, 0.67, 1.13, 2.10, 0.44,
        1.70, 0.61, 0.43, 0.48, 0.47, 0.31])
    #ax.errorbar(dt, f*scale_sn2003bg, ef*scale_sn2003bg, fmt='o', c=cols[0],
    #        label)
    ax.plot(dt, f*scale_sn2003bg, c=cols[0], label='03bg')

    # 15 GHz
    ax = axarr[1,0]
    ax.set_title("15 GHz")

    # SN2003L
    dt = np.array([33.3, 41.6, 44.5, 46.4, 52.3, 55.3, 56.3, 60.3, 66.4, 68.3,
        72.4, 75.4, 81.3, 99.4, 135.3, 147, 167.1, 189, 204, 226.8, 257.7,
        283.7, 299.7, 328.5, 339.6, 476])
    f = (1/1000)*np.array([2287, 2686, 2753, 2719, 3014, 2799, 2865, 2924, 2857,
        2966, 2901, 3055, 2664, 2581, 1925, 1452, 1279, 1135, 1177, 824, 768,
        612, 520, 148, 506, 125])
    ef = (1/1000)*np.array([121, 156, 153, 144, 137, 149, 131, 132, 123, 115,
        150, 173, 118, 320, 282, 235, 226, 125, 166, 237, 165, 122, 124, 148,
        149, 137])
    #ax.errorbar(dt, f*scale_sn2003L, ef*scale_sn2003L, fmt='o', c=cols[1])
    ax.plot(dt, f*scale_sn2003L, c=cols[1], label='03L', ls='--')

    # AT2020xnd
    dt = np.array([18, 25, 71.9, 94.7, 131.6])
    f = np.array([0.037, 0.095, 0.401, 0.278, 0.122])
    ef = np.array([0.01, 0.011, 0.046, 0.032, 0.016])
    ax.errorbar(dt, f*scale_at2020xnd, ef*scale_at2020xnd, fmt='o', c=cols[2],
            zorder=20, label='20xnd')

    # SN 2003bg
    dt = np.array([23, 35, 48, 58, 63, 73, 85, 91, 115,
        129, 132, 142, 157, 161, 181, 201, 214, 227, 255])
    f = np.array([47, 62.11, 71.51, 69.31, 64.75, 41.18, 19.39, 17.34,
        32.14, 41.02, 40.92, 39.24, 30.88, 29.17, 21.06, 16.58, 16.54, 15.58,
        12.72])
    ef = np.array([0.87, 1.26, 1.46, 1.41, 1.32, 0.87, 0.54, 0.83, 0.68, 0.84,
        0.83, 0.81, 0.63, 0.61, 0.48, 0.44, 0.39, 0.35, 0.35])
    #ax.errorbar(dt, f*scale_sn2003bg, ef*scale_sn2003bg, fmt='o', c=cols[0])
    ax.plot(dt, f*scale_sn2003bg, c=cols[0], label='03bg')

    # 10 GHz (or 8.5)
    ax = axarr[1,1]
    ax.set_title("10 GHz")

    # AT2020xnd
    dt = np.array([13, 25, 37, 52, 72, 95, 131.6])
    f = np.array([0.024, 0.057, 0.079, 0.154, 0.180, 0.168, 0.109])
    ef = np.array([0.006, 0.005, 0.010, 0.005, 0.023, 0.022, 0.01])
    ax.errorbar(dt, f*scale_at2020xnd, ef*scale_at2020xnd, fmt='o', c=cols[2],
            zorder=20, label='20xnd')

    # SN2003bg
    dt = np.array([10, 12, 23, 35, 48, 58, 63, 73, 85, 91, 115,
        129, 132, 142, 157, 161, 181, 201, 214, 227, 242])
    f = np.array([2.51, 3.86, 12.19, 24.72, 40.34, 51.72, 49.64,
        46.20, 38.65, 33.85, 45.74, 53.94, 54.27, 54.83, 48.43,
        47.43, 35.76, 31.35, 28.67, 27.38, 24.57])
    ef = np.array([0.07, 0.1, 0.26, -.5, 0.81, 1.04, 1, 0.93, 0.79,
        0.71, 0.92, 1.08, 1.09, 1.10, 0.97, 0.95, 0.72, 0.63, 0.58,
        0.55, 0.50])
    #ax.errorbar(dt, f*scale_sn2003bg, ef*scale_sn2003bg, fmt='o', c=cols[0])
    ax.plot(dt, f*scale_sn2003bg, c=cols[0], label='03bg')

    # SN2003L
    dt = np.array([25,27.3,28.3,29.4,30.3,31.2,32.3,36.3,38.3,41.6,44.5,46.4,
        48.4,53.3,55.3,60.3,66.4,68.3,72.4,75.4,78.3,81.3,85.4,90.2,119.0])
    f = (1/1000)*np.array([743,810,848,966,1051,947,883,1147,1211,1483,
        1448,1413,1546,1854,1969,2283,2351,2081,2527,2673,2422,2500,2776,2551,2532])
    ef = (1/1000)*np.array([39,63,65,60,62,52,53,47,39,76,62,55,45,41,53,51,47,
        52,52,57,61,43,51,54,75])
    #print(len(dt), len(f), len(ef))
    #ax.errorbar(dt, f*scale_sn2003L, ef*scale_sn2003L, fmt='o', c=cols[1])
    ax.plot(dt, f*scale_sn2003L, c=cols[1], label='03L', ls='--')

    # Formatting
    for ax in axarr.flatten():
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(10,140)
        ax.set_ylim(0.006,0.8)
        ax.legend(loc='lower right', ncol=2, fontsize=8)

    axarr[0,0].set_ylabel("Flux Density (mJy)")
    axarr[1,0].set_ylabel("Flux Density (mJy)")
    axarr[1,0].set_xlabel("Days")
    axarr[1,1].set_xlabel("Days")

    plt.tight_layout()
    #plt.show()
    plt.savefig("compare_lc.png", dpi=300)
    plt.close()


def sed():
    """ Look at the SEDs and measure the power-law index """
    fig, axarr = plt.subplots(2,2,figsize=(6,6))
    d = get_format()
    cols = d['colors']['4']

    # SN2003bg: dt=23d
    ax = axarr[0,0]
    ax.set_title("SN03bg, 23d")

    nu = np.array([1.43, 4.86, 8.46, 15, 22.5])
    f = np.array([0.55, 2.75, 12.19, 47.03, 106.30])
    ef = np.array([0.20, 0.11, 0.26, 0.98, 2.16])
    ax.errorbar(nu[1:], f[1:], ef[1:], fmt='o', c=cols[0])
    ax.scatter(nu[0], f[0], marker='v', c=cols[0])

    xplt = np.linspace(1,50)
    yplt = f[1]*(xplt/nu[1])**(2.5)
    ax.plot(xplt, yplt, c=cols[0], lw=0.5, label='$\\nu^{5/2}$')
    yplt = f[1]*(xplt/nu[1])**(2)
    ax.plot(xplt, yplt, c=cols[0], lw=0.5, ls='--', label='$\\nu^2$')
    ax.set_xlim(1.1, 30)
    ax.set_ylim(0.3, 170)
    ax.legend()

    # SN2003L: 
    ax = axarr[0,1]
    ax.set_title("SN03L, 33d")

    nu = np.array([8.5, 15, 22.5])
    f = np.array([883, 2287, 2973])/1000
    ef = np.array([53, 121, 99])/1000
    ax.errorbar(nu, f, ef, fmt='o', c=cols[1])

    xplt = np.linspace(1,50)
    yplt = f[1]*(xplt/nu[1])**(2.5)
    ax.plot(xplt, yplt, c=cols[1], lw=0.5, label='$\\nu^{5/2}$')
    yplt = f[1]*(xplt/nu[1])**(2)
    ax.plot(xplt, yplt, c=cols[1], lw=0.5, ls='--', label='$\\nu^2$')
    ax.set_xlim(7, 26)
    ax.set_ylim(0.7, 3.6)
    ax.legend()

    # SN2003L late-time
    ax = axarr[1,1]
    ax.set_title("SN03L, 147d")

    nu = np.array([4.9, 8.5, 15, 22.5])
    f = np.array([2380, 2048, 1177, 670])/1000
    ef = np.array([81, 66, 166, 148])/1000
    ax.errorbar(nu, f, ef, fmt='o', c=cols[1])

    xplt = np.linspace(1,50)
    yplt = f[1]*(xplt/nu[1])**(-1)
    ax.plot(xplt, yplt, c=cols[1], lw=0.5, label='$\\nu^{-1}$')
    yplt = f[1]*(xplt/nu[1])**(-2)
    ax.plot(xplt, yplt, c=cols[1], lw=0.5, ls='--', label='$\\nu^{-2}$')
    ax.set_xlim(3, 26)
    ax.set_ylim(0.4, 3)
    ax.legend()

    # SN2003bg late-time
    ax = axarr[1,0]
    ax.set_title("SN03bg, 58d")
    nu = np.array([4.86, 8.46, 15, 22.5, 43.3])
    f = np.array([22.37, 51.72, 69.31, 39.91, 20.56])
    ef = np.array([0.46, 1.04, 1.41, 0.89, 5.87])
    ax.errorbar(nu, f, ef, fmt='o', c=cols[0])

    xplt = np.linspace(1,50)
    yplt = f[3]*(xplt/nu[3])**(-1)
    ax.plot(xplt, yplt, c=cols[0], lw=0.5, label='$\\nu^{-1}$')
    yplt = f[3]*(xplt/nu[3])**(-2)
    ax.plot(xplt, yplt, c=cols[0], lw=0.5, ls='--', label='$\\nu^{-2}$')
    ax.set_xlim(4, 56)
    ax.set_ylim(13, 82)
    ax.legend()

    for ax in axarr.flatten():
        ax.set_yscale('log')
        ax.set_xscale('log')

    axarr[0,0].set_ylabel("Flux Density (mJy)")
    axarr[1,0].set_ylabel("Flux Density (mJy)")
    axarr[1,0].set_xlabel("Freq [GHz]")
    axarr[1,1].set_xlabel("Freq [GHz]")

    plt.tight_layout()
    #plt.show()
    plt.savefig("sed.png", dpi=300)
    plt.close()


if __name__=="__main__":
    sed()
