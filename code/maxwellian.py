""" Fit a relativistic Maxwellian """

import matplotlib.pyplot as plt
from get_radio import *
from scipy.optimize import curve_fit


islim, tel, freq, days, flux, eflux = get_data_all()
bin = 71
col = 'cornflowerblue'

fig,ax = plt.subplots(1,1,figsize=(6,4))

choose = np.logical_and.reduce((
    days>bin-bin/20, days<bin+bin/20, eflux>0))
if bin!=71:
    order = np.argsort(freq[choose])
    ax.errorbar(
            freq[choose][order], flux[choose][order], eflux[choose][order],
            fmt='%s-' %'D', c=col, label='%s d' %bin)
    nondet = np.logical_and(choose, eflux==0)
    #for i,x in enumerate(freq[nondet]):
    #    print(x)
    #    y = flux[nondet][i]
    #    print(y)
    #    ax.arrow(x, y,
    #            0, -y/2, head_width=x/7,
    #            length_includes_head=True, head_length=y/7, color=col)

else:
    keep = freq[choose] < 80
    order = np.argsort(freq[choose][keep])
    ax.errorbar(
            freq[choose][keep][order], flux[choose][keep][order], 
            eflux[choose][keep],
            fmt='%s-' %'D', c=col, label='%s d' %bin)


def fitfunc(nu_ghz, A, B, T):
    nu = nu_ghz * 1E9
    m_e = 9.1E-28
    c = 3E10
    e = 4.8E-10
    k = 1.38E-16

    theta_e = k*T / (m_e * c**2)
    w = nu*(2*np.pi)
    w_b = e*B/(m_e*c)
    X = w/w_b
    xM = 2*X/(3*theta_e**2)

    y = A*2.5651*(1+1.92/xM**(1/3)+0.9977/xM**(2/3))*np.exp(-1.8899*xM**(1/3))
    return y


p0 = [1000, 10, 3E10]
popt, pcov = curve_fit(
        fitfunc, freq[choose][order][1:], flux[choose][order][1:], 
        p0, sigma=eflux[choose][order][1:], absolute_sigma=True)
# Best fit: A = 23.9 +/- 19.2
# B and T basically unconstrained...huge uncertainties
# B = 19.64 G
# T = 4.32E10 K

nufit = np.arange(50, 200) 

yfit = fitfunc(nufit, *popt)
ax.plot(nufit,yfit, c='k', ls='--', zorder=5)

nufit = np.arange(20, 100)
yfit = flux[choose][order][1]*((nufit)/freq[choose][order][1])**(1/3)
ax.plot(nufit,yfit, ls=':', c='grey')
ax.text(30, 0.7, r'$\nu^{1/3}$')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel("Frequency [GHz]")
ax.set_ylabel("Flux Density (mJy)")

ax.text(100, 0.8, r'$\propto e^{-1.8899\,x_M^{1/3}}$', fontsize=12)


#ax.set_ylim(1E-1, 1)

plt.tight_layout()
plt.show()
#plt.savefig("camel_sed_maxwellian.png", dpi=300)
