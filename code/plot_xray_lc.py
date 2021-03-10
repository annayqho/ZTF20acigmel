""" Plot the X-ray light curve compared to that of AT2018cow """
import matplotlib
from astropy.time import Time
from get_radio import *

# Font sizes
large = 14
medium = 11
small = 9

# Redshift for transforming to the rest-frame
z = 0.2442

fig,ax = plt.subplots(1,1,figsize=(6,4))

dt = np.array([25.8, 31.8, 47.1, 75.2])/(1+z)
ymin = np.array([1.27, 0.67, 0.11])*10
ymax = np.array([0.96, 0.75, 0.17])*10
y = np.array([3.46, 2.79, 0.15, 0.22])*10
ax.errorbar(dt[0:3], y[0:3], yerr=(ymin, ymax), c='k', fmt='o',
        label='AT2020xnd')

ax.scatter(dt[-1], y[-1], marker='_', c='k', s=50)
ax.arrow(dt[-1], y[-1], 0, -y[-1]/3, color='k',
        length_includes_head=True, head_length=y[-1]/7,
        head_width=dt[-1]/30)

# Overplot the 18cow X-ray light curve
t0_offset = (Time('2018-06-19T10:34:30.742') - Time(58285, format='mjd')).value
dat = Table.read("../data/AT2018cow/swift_lc.txt", format='ascii')
t = dat['col1'].data / (3600*24) # in days
# count to flux conversion (absorbed): 4.26 × 10-11 erg cm-2 ct-1
flux = dat['col4'].data * 4.26 * 10**(-11)
eflux = dat['col5'].data * 4.26 * 10**(-11)
t_xray = (t-t[0])+3

# Overplot the CSS161010 X-ray data
ax.errorbar([99, 130], [1.33E-15, 1.94E-15], [0.76E-15, 0.97E-15],
        fmt='s-', label='CSS161010')

ax.plot(t_xray, flux/1E-15/400, lw=0.5, c='k', label='AT2018cow')

ax.set_xlabel("Rest-frame days since $t_0$", fontsize=large)
ax.set_ylabel("$10^{-15}$ erg cm$^{-2}$ s$^{-1}$", fontsize=large)
ax.set_ylim(2E-1,65)
ax.set_yscale('log')
ax.set_yticks([0.2, 1, 10])
ax.set_xticks([10, 20, 30, 50, 70])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.legend(loc='upper right', fontsize=medium)
ax.minorticks_off()
ax.tick_params(axis='both', labelsize=large)

plt.tight_layout()
plt.show()
#plt.savefig("xray_lc.png", dpi=200)
#plt.close()
