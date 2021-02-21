""" Plot radio/mm light curve """

import matplotlib
from astropy.time import Time
from get_radio import *

# Font sizes
large = 14
medium = 11
small = 9

# Marker sizes
ms=10

# Redshift for transforming to the rest-frame
z = 0.2442

# Initializing plot
fig,axarr = plt.subplots(2,1,figsize=(6,8),gridspec_kw={'height_ratios': [2,1]}, sharex=True)

def plot_radio(ax):
    # Get radio data
    islim, tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()

    # Plot the NOEMA LC
    choose = freq==79
    ax.errorbar(
            days[choose]/(1+z), flux[choose], 
            yerr=np.sqrt(eflux_form[choose]**2+eflux_sys[choose]**2), 
            fmt='o', c='k', label='NOEMA 100 GHz', ms=10)

    # Plot line of t^2
    xvals = np.linspace(10,40)
    yvals = flux[choose][1]*(xvals/days[choose][1])**2
    ax.plot(xvals/(1+z),yvals,c='k',lw=0.5,ls='--')
    ax.text(17, 0.5, '$f_\\nu \propto t^2$', fontsize=large,
            ha='right')

    # Plot the ATCA 34 GHz LC
    # also use VLA 33 GHz
    choose = np.logical_and(
            np.logical_or(freq==27.3,freq==26.5), islim==False)
    ax.errorbar(
            days[choose]/(1+z), flux[choose], 
            yerr=np.sqrt(eflux_form[choose]**2+eflux_sys[choose]**2), 
            fmt='D', c='orange', label='ATCA/VLA 34 GHz', ms=10)

    choose = np.logical_and(freq==27.3, islim==True)
    ax.scatter(days[choose][0]/(1+z), flux[choose][0], marker='_', c='orange', s=ms*10)
    ax.arrow(
            days[choose][0]/(1+z), flux[choose][0], 0, -flux[choose][0]/5, 
            length_includes_head=True, head_length=flux[choose][0]/10, 
            head_width=days[choose][0]/20, color='orange')

    # Plot the VLA 10 GHz LC
    choose = np.logical_and(freq==8, islim==False)
    ax.errorbar(
            days[choose]/(1+z), flux[choose], 
            yerr=np.sqrt(eflux_form[choose]**2+eflux_sys[choose]**2), 
            fmt='s', c='purple', label='VLA 10 GHz', ms=10)


def plot_xray(ax):
    dt = np.array([26, 32, 48])/(1+z)
    ymin = np.array([0.93,0.71,0.01])
    ymax = np.array([1.94,1.61,0.55])
    y = np.array([1.37, 1.10, 0.26])
    ax.errorbar(dt, y, yerr=(y-ymin, ymax-y), c='k', fmt='o')

    # Overplot the 18cow X-ray light curve
    t0_offset = (Time('2018-06-19T10:34:30.742') - Time(58285, format='mjd')).value
    dat = Table.read("../data/AT2018cow/swift_lc.txt", format='ascii')
    t = dat['col1'].data / (3600*24) # in days
    # count to flux conversion (absorbed): 4.26 × 10-11 erg cm-2 ct-1
    flux = dat['col4'].data * 4.26 * 10**(-11)
    eflux = dat['col5'].data * 4.26 * 10**(-11)
    t_xray = (t-t[0])+3

    ax.plot(t_xray, flux/1E-15/400, lw=0.5, c='k', label='AT2018cow')


plot_radio(axarr[0])
plot_xray(axarr[1])

# Formatting
ax = axarr[0]
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=large)
ax.set_ylim(0.015,1.2)
ax.set_xlim(11,150)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_yticks([0.1, 1, 0.05])
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.legend(loc='lower right', fontsize=medium)
ax.minorticks_off()

#ax.axvline(x=124/1.2442)

ax = axarr[1]
ax.set_xlabel("Rest-frame days since $t_0$", fontsize=large)
ax.set_ylabel("$10^{-15}$ erg cm$^{-2}$ s$^{-1}$", fontsize=large)
ax.set_ylim(1E-1,40)
ax.set_yscale('log')
ax.set_yticks([0.1, 1, 10])
ax.set_xticks([15, 20, 30, 50, 70])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.legend(loc='upper right', fontsize=medium)
ax.minorticks_off()

for ax in axarr:
    ax.tick_params(axis='both', labelsize=large)

# Display
plt.tight_layout()
#plt.show()
plt.savefig("radio_lc.png", dpi=300)
