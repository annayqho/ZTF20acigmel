""" Plot radio/mm light curves
If you plot frequencies with >=3 detections then you have
8, 12, 17.7, 26.5, 36.2, 79, 94 (7 bands)

If you plot frequencies with >=2 detections then you add
4, 14.5, 27.3, 131 (4 bands; total 11)

I think you should plot frequencies with >=3 OBSERVATIONS,
because those are important too.
"""

import matplotlib
from astropy.time import Time
from get_radio import *
from ssa_lc_fixalpha import *
from format import *

form = get_format()
large = form['font_large']
medium = form['font_med']
small = form['font_small']

# Marker sizes
ms=10

# Initializing plot
fig,axarr = plt.subplots(4,3,figsize=(7.5,8), sharex=True, sharey=True)


def plot_all(ax):
    """ Plot all of the frequencies as grey in the background """
    islim, tel, freq, days, flux, eflux = get_data_all()
    plot_freqs = np.unique(freq)
    plot_freqs = plot_freqs[plot_freqs < 240]

    for plot_freq in plot_freqs:
        choose = np.logical_and(freq==plot_freq, islim==False)
        ax.plot(days[choose], flux[choose], c='grey', alpha=0.3)

    # Plot a cross marking the peak of the 79 GHz light curve
    ax.scatter(31.8, 1.076, marker='+', c=form['colors']['3'][1], zorder=10)

    # Plot vertical lines indicating epochs for SEDs
    for b in [18, 24, 30.3, 38, 46, 51.9, 71, 95, 132]:
        ax.axvline(x=b, lw=0.5, ls=':')


def plot_panel(ax, choose):
    """ Plot a single panel """
    islim, tel, freq, days, flux, eflux = get_data_all()

    # Plot LC at a single frequency
    det = np.logical_and(islim==False, choose)
    ax.errorbar(
            days[det], flux[det], yerr=eflux[det],
            fmt='o', c='k', ms=5)
    nondet = np.logical_and(islim==True, choose)
    ax.errorbar(
            days[nondet], flux[nondet], yerr=eflux[nondet],
            fmt='o', c='k', mfc='white', mec='k', ms=5)

    ax.plot(days[choose], flux[choose], c='k')

    for ii,x in enumerate(days[nondet]):
        y = flux[nondet][ii]
        ax.arrow(x, y, 0, -y/3, head_width=x/10,
            length_includes_head=True, head_length=y/10, color='k')


def plot_tpow(ax, choose, pow, ind=1, tjust='right'):
    """ Plot line of t^2 """
    islim, tel, freq, days, flux, eflux = get_data_all()

    xvals = np.linspace(10,200)
    yvals = flux[choose][ind]*(xvals/days[choose][ind])**pow
    ax.plot(xvals,yvals,c='k',lw=0.5,ls='--')
    ax.text(days[choose][ind], flux[choose][ind], 
            '$\propto t^{%s}$' %pow, 
            fontsize=small, ha=tjust, va='bottom')


 
if __name__=="__main__":
    islim, tel, freq, days, flux, eflux = get_data_all()

    for ax in axarr.flatten():
        plot_all(ax)

    plot_freq = []

    # Plot all frequencies with >=2 observations
    for f in np.unique(freq):
        print(f)
        if np.logical_and(sum(freq==f)>=2, f!=34):
            plot_freq.append(f)
    plot_freq = np.array(plot_freq)[::-1]
    print(plot_freq)

    for ii,f in enumerate(plot_freq):
        ax = axarr.flatten()[ii]
        choose = freq==f
        fstr = int(f)
        if f==33:
            choose = np.logical_and(freq>32, freq<35)
            fstr = "33/34"
        if f==230:
            choose = np.logical_and(freq>225, freq<235)
            fstr = "227/230"
        plot_panel(ax, choose)
        ax.text(
                0.95,0.15,'%s GHz' %fstr,fontsize=medium,
                transform=ax.transAxes,ha='right', va='top')

        if np.logical_or.reduce((f==117, f==98, f==33, f==15)):
            plot_tpow(ax, choose, 2)
        if f==10:
            plot_tpow(ax, choose, 1, ind=2)
        if f==94.244:
            plot_tpow(ax, choose, 2, 1, 'right')
            plot_tpow(ax, choose, -4, 5, 'left')
        if f==78.756:
            plot_tpow(ax, choose, 2, 1, 'right')
            plot_tpow(ax, choose, -4, 5, 'left')
        if f==45:
            plot_tpow(ax, choose, -4, 2, 'left')
        if f==33:
            plot_tpow(ax, choose, -4, 6, 'left')
        if f==22:
            plot_tpow(ax, choose, -4, 1, 'left')
        if f==18:
            plot_tpow(ax, choose, 1, 1, 'right')

    # Formatting
    axarr[0,0].set_yscale('log')
    axarr[0,0].set_ylim(0.015,1.3)
    axarr[0,0].set_xscale('log')
    axarr[0,0].set_xlim(9,170)
    for ax in axarr[:,0]:
        ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=large)
        ax.tick_params(axis='both', labelsize=large)
    for ax in axarr[-1,:]:
        ax.set_xlabel("$\Delta t$ [d]", fontsize=large)
        ax.tick_params(axis='both', labelsize=large)
    #axarr[1,2].set_xlabel("$\Delta t$ [d]", fontsize=large)
    #axarr[1,2].tick_params(axis='both', labelsize=large)
    #axarr[1,3].set_xlabel("$\Delta t$ [d]", fontsize=large)
    #axarr[1,3].tick_params(axis='both', labelsize=large)

    # Display
    #plt.tight_layout()
    plt.subplots_adjust(hspace=0.05, wspace=0.1)
    #plt.show()
    plt.savefig("radio_lc_with_epochs.png", dpi=300, bbox_inches='tight')
    plt.close()
