""" Fit a relativistic Maxwellian """

import matplotlib.pyplot as plt
from get_radio import *
from plot_radio_sed import day46

day46(days, freq, flux, eflux, ax):

fig,ax = plt.subplots(1,1,figsize=(6,4))

islim, tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
day46(days, freq, flux, eflux, ax)

m_e = 9.1E-28
c = 3E10
e = 4.8E-10
k = 1.38E-16

B = 10
T = 3E10

nufit = np.arange(68, 178) * 1E9

theta_e = k*T / (m_e * c**2)
w = nufit*(2*np.pi)
w_b = e*B/(m_e*c)
X = w/w_b
xM = 2*X/(3*theta_e**2)

yfit = 2.5651*np.exp(-1.8899*xM**(1/3))
ax.plot(nufit/1E9,1000*yfit)

nufit = np.arange(20, 100)
yfit = 0.8*(nufit/78)**(1/3)
ax.plot(nufit,yfit)

nufit = np.arange(8, 30)
yfit = 0.13*(nufit/8)**(2)
ax.plot(nufit,yfit)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel("Frequency [GHz]")
ax.set_ylabel("Flux Density (mJy)")

ax.text(100, 0.8, r'$\propto e^{-1.8899\,x_M^{1/3}}$', fontsize=12)

#ax.set_ylim(1E-1, 1)

plt.tight_layout()
plt.show()
#plt.savefig("camel_sed_maxwellian.png", dpi=300)
