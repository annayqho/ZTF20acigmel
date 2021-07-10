""" Fit a relativistic Maxwellian """

import matplotlib.pyplot as plt
from get_radio import *
from scipy.optimize import curve_fit


def fitfunc(nu_ghz, A, C):
    nu = nu_ghz * 1E9
    m_e = 9.1E-28
    c = 3E10
    e = 4.8E-10
    k = 1.38E-16

    #theta_e = k*T / (m_e * c**2)
    #w_b = e*B/(m_e*c)
    #X = w/w_b
    #xM = 2*X/(3*theta_e**2)

    w = nu*(2*np.pi)
    xM = C*(2*m_e**3*c**5*w) / (3*k**2*e)
    y = A*2.5651*(1+1.92/xM**(1/3)+0.9977/xM**(2/3))*np.exp(-1.8899*xM**(1/3))

    return y


islim, tel, freq, days, flux, eflux = get_data_all()
bin = 46
col = '#a05195'
marker = '*'
msize = 14

fig,ax = plt.subplots(1,1,figsize=(5.5,4))

choose = np.logical_and.reduce((
    days>bin-bin/20, days<bin+bin/20, eflux>0))
if bin!=71:
    order = np.argsort(freq[choose])
    ax.errorbar(
            freq[choose][order], flux[choose][order], eflux[choose][order],
            fmt='%s-' %marker, c=col, label=None, ms=msize)
    ax.scatter(0, 0, marker='%s' %marker, c=col, label='%s d' %bin, s=100)
    nondet = np.logical_and(choose, eflux==0)
else:
    keep = freq[choose] < 80
    order = np.argsort(freq[choose][keep])
    ax.errorbar(
            freq[choose][keep][order], flux[choose][keep][order], 
            eflux[choose][keep],
            fmt='%s-' %'D', c=col, label='%s d' %bin)


p0 = [23.9, 3E-23]
popt, pcov = curve_fit(
        fitfunc, freq[choose][order][1:], flux[choose][order][1:], 
        p0, sigma=eflux[choose][order][1:], absolute_sigma=True)

# I actually get (2.73 +/- 1.13) x 10^(-23)

B = np.logspace(0,2)
T = np.sqrt(3.7E22 / B)

nu_g = e*B / (m_e*c)
nu_cutoff = (k*T/(m_e*c**2))**2 * nu_g

# Best fit: A = 23.9 +/- 19.2
# B and T basically unconstrained...huge uncertainties
# oh, that's because there was a degeneracy between them
# B = 19.64 G
# T = 4.32E10 K

nufit = np.arange(50, 300) 

yfit = fitfunc(nufit, *p0)
ax.plot(nufit,yfit, c='k', ls='--', zorder=5, label='Maxwellian')

# Also plot a power law
yfit = flux[choose][order][1]*(nufit/freq[choose][order][1])**(-1.5)
ax.plot(nufit,yfit, c='k', ls=':', zorder=5, label=r'Power law ($\nu^{-1.5}$)')

yfit = flux[choose][order][1]*(nufit/freq[choose][order][1])**(-1)
ax.plot(nufit,yfit, c='k', ls='-', lw=0.8, zorder=5, label=r'Power law ($\nu^{-1}$)')

print(freq[choose][order], flux[choose][order])


ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xticks([30,50,100,200])
ax.set_xticklabels([30,50,100,200])
ax.set_yticks([0.1, 0.2, 0.3, 0.4, 0.6, 1])
ax.set_yticklabels([0.1, 0.2, 0.3, 0.4, 0.6, 1])
plt.minorticks_off()

ax.set_xlabel("Frequency [GHz]", fontsize=14)
ax.set_ylabel("Flux Density (mJy)", fontsize=14)
ax.tick_params(axis='both', labelsize=12)

ax.set_xlim(30, 250)
ax.set_ylim(1E-1, 1)

plt.legend(loc='lower left', fontsize=11)
plt.tight_layout()
#plt.show()
plt.savefig("camel_sed_maxwellian.png", dpi=300)
plt.close()
