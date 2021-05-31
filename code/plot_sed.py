""" Plot full radio to X-ray SED at dt = 30 +/- 3 days """

import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15


bigsize=16
smallsize=10

def plot_xray(flux, c):
    """ Take the integrated flux across 0.3-10 keV,
    and use that and the geometric mean of the frequency and the spectral index
    to solve for the normalization coefficient for your spectrum

    Parameters
    ----------
    flux: flux
    c: color of the point/line

    2.4E18 Hz corresponds to 10 keV
    7.2E16 Hz corresponds to 0.3 keV
    """
    nu0 = np.sqrt(0.3*10) * (2.4E18/10)
    alpha = 0.75 # photon index is 1.75 at 32d
    nu2 = 2.4E18
    nu1 = 7.2E16
    A = flux*(1-alpha) / (nu0**alpha * (nu2**(1-alpha) - nu1**(1-alpha)))
    xplot = np.linspace(7.2E16, 2.4E18, 1000)
    print(A)
    print(nu0)
    yplot = A*(xplot/nu0)**(-alpha)
    plt.plot(xplot, xplot*yplot, c=c, ls='-', label='Chandra')


fig,ax = plt.subplots(1,1,figsize=(5,3.5), dpi=200)
dcm = Planck15.luminosity_distance(z=0.2442).cgs.value

# from Perley+2021
# at dt=26d ... MJD 59158
opt_nu = 3E18 / np.array([3904, 5055, 6193]) # in Hz
mAB = np.array([23.56, 23.44, 23.52])
fnu = 10**((mAB+48.6)/(-2.5))
opt_Lnu = 4*np.pi*dcm**2 * fnu
plt.scatter(opt_nu, opt_nu*opt_Lnu, c='k', marker='s', label='VLT+FORS2')

# this is 0.3-10 keV (from Yuhan's analysis)
# on MJD 59158
xray_L = 3.46E-14*4*np.pi*dcm**2
plot_xray(xray_L, 'k')

# In the rest-frame
radio_nu = np.array([79,94,12,8,4])*1E9
radio_L = np.array([0.679,0.648,0.095,0.057,0.046])*1E-3*1E-23*4*np.pi*dcm**2
plt.scatter(radio_nu, radio_nu*radio_L, c='k', marker='D', label='VLA/NOEMA')

# Plot the radio to X-ray spectral index
xplot = np.linspace(1E11, 1E17)
yplot = (1.1E41)*(xplot/1E11)**0.18
plt.plot(xplot, yplot, lw=0.5, ls='--', c='k')
plt.text(1E14, 2E41, r'$f_\nu\propto \nu^{0.18}$', fontsize=smallsize)

plt.text(0.05, 0.95, '20d after explosion (rest-frame)', fontsize=smallsize,
        transform=ax.transAxes)

plt.xlim(1E9, 5E18)
plt.ylim(1E38, 1E43)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Frequency [Hz]", fontsize=bigsize)
plt.ylabel(r"$\nu L_\nu$ (erg/s)", fontsize=bigsize)
plt.xticks(fontsize=bigsize)
plt.yticks(fontsize=bigsize)
plt.legend(fontsize=smallsize, loc='lower right')

plt.tight_layout()
#plt.savefig("sed.png")
#plt.close()

plt.show()
