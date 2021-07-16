""" Fit a relativistic Maxwellian """

import matplotlib.pyplot as plt
from get_radio import *
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15


def fitfunc(nu_ghz, L, tau_m, nu_m):
    """ Treat all frequencies in GHz """
    # Eq 29 from Mahadevan paper
    xM = nu_ghz/nu_m
    Inu = 2.5651*(1+1.92/xM**(1/3)+0.9977/xM**(2/3))*np.exp(-1.8899*xM**(1/3))

    # Ben's equation
    exponent = -tau_m * (nu_ghz/nu_m)**(-1) * Inu
    Lnu = L * xM**2 * (1-np.exp(exponent))

    return Lnu


def fitfunc_physical(nu_ghz, B, ne):
    """ This time with physical parameters

    Te_rel = kT/mc^2
    beta = v/c
    R = in cm
    B = in Gauss

    Treat all frequencies in GHz """
    beta = 0.2

    # Scale factor in mJy
    v_cgs = beta*3E10
    L = 4.59E22 * B**2 * (beta/0.1)**12 * (td/50)**2
    dcm = Planck15.luminosity_distance(z=0.2442).cgs.value
    fmjy = (L / (4*np.pi*dcm**2)) / 1E-23 / 1E-3

    # Synchrotron frequency in GHz
    nu_m = (0.033) * (beta/0.1)**4 * B

    # taum
    tau_m = 1.18E6 * ne * B**(-1) * (beta/0.1)**(-9) * (td/50)

    # Eq 29 from Mahadevan paper
    xM = nu_ghz/nu_m
    Inu = 2.5651*(1+1.92/xM**(1/3)+0.9977/xM**(2/3))*np.exp(-1.8899*xM**(1/3))

    # Ben's equation
    exponent = -tau_m * (nu_ghz/nu_m)**(-1) * Inu
    fnu = fmjy * xM**2 * (1-np.exp(exponent))

    return fnu


def fitfunc_powlaw(x, A, beta):
    return A*x**(beta)

fig,ax = plt.subplots(1,1,figsize=(5.5,4))

islim, tel, freq, days, flux, eflux = get_data_all()
z = 0.2442

# Choose two epochs of observations to fit together
#choose = np.logical_and.reduce((days>41, days<52, islim==False))
choose = np.logical_and.reduce((days>70, days<72, islim==False))
# Get values in rest-frame
x = freq[choose] * (1+z)
y = flux[choose] / (1+z)
ey = eflux[choose] / (1+z)
# Sort
order = np.argsort(x)
x = x[order]
y = y[order]
ey = ey[order]
td = np.average(days[choose])

# Plot the data
#col = '#a05195'
col = 'darkblue'
#marker = '*'
marker = 'o'
#msize = 14
msize = 10
ax.errorbar(x, y, ey, 
        fmt='%s-' %marker, c=col, label=None, ms=msize)

# Fit for a Maxwellian w/o physical parameters
p0 = [5E-4, 5E4, 1]

# Fit for the Maxwellian in terms of physical quantities
p0 = [2.0, 40.0]
popt, pcov = curve_fit(fitfunc_physical, x[1:], y[1:], sigma=ey[1:], absolute_sigma=True, p0=p0)
xfit = np.linspace(1,300)
yfit = fitfunc_physical(xfit, *popt)
print("Maxwellian fit:")
for i,param in enumerate(popt):
    print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
ax.plot(xfit,yfit, c='k', ls='--', zorder=5, label='Maxwellian')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xticks([10,30,50,100,200])
ax.set_xticklabels([10,30,50,100,200])
ax.set_yticks([0.1, 0.2, 0.3, 0.4, 0.6, 1])
ax.set_yticklabels([0.1, 0.2, 0.3, 0.4, 0.6, 1])
plt.minorticks_off()

ax.set_xlabel("Rest Frequency [GHz]", fontsize=14)
ax.set_ylabel("Rest Flux Density (mJy)", fontsize=14)
ax.tick_params(axis='both', labelsize=12)

ax.set_xlim(9, 300)
ax.set_ylim(7E-2, 1)

plt.legend(loc='upper left', fontsize=11)
plt.tight_layout()
plt.show()
#plt.savefig("camel_sed_maxwellian.png", dpi=300)
#plt.close()
