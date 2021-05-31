""" Evolution of physical parameters, using the SSA model """
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from astropy.cosmology import Planck15
from ssa_lc import *

eps = 1
eps_B = 1/3
f = 0.5

D_mpc = Planck15.luminosity_distance(z=0.2442).value
t0 = 28
Fa0 = 2.06
nua0 = 95.6
alpha_r=1.54
k = 2.63

m_p = 1.67E-24

t = np.logspace(1,2,1000)
t_sec = t*86400
out = Fnu(10, t, alpha_r, k, s=1.00, opt_thick_index=1.12, t0=t0, Fa0=Fa0, nua0=nua0)
F = out[0]
is_right_1 = out[1]

B0 = 0.30 * eps**(-4/19) * (Fa0/1000)**(-2/19) * (D_mpc)**(-4/19) * (nua0/5)
R0 = 5.1E15*eps**(-1/19)*(Fa0/1000)**(9/19)*(D_mpc)**(18/19)*(nua0/5)**(-1)
B = B0*(t/t0)**(-alpha_r*k/2)
R = R0*(t/t0)**(alpha_r)
v = R/t_sec
E = (1/eps_B)*(4*np.pi/3)*f*R**3*B**2/(8*np.pi)
ne = (4/3) * (1/m_p) * (1/eps_B) * (B**2/(8*np.pi)) * (1/v**2)
vw = 1000E5
Mdot = 4*np.pi*m_p * ne * R**2 * vw # grams per second
Mdot = Mdot / 2E33 # solar masses per second
Mdot = Mdot * 86400 * 365 # solar masses per year

fig,axarr = plt.subplots(3,2,figsize=(5,6))

ax = axarr[0,0]
ax.plot(t,E/1E48,c='k')
ax.set_ylabel("Energy ($10^{48}$ erg)")
ax.text(0.05, 0.95, "$E \propto t^{0.56}$", fontsize=12, transform=ax.transAxes,
        ha='left', va='top')

ax = axarr[1,0]
ax.plot(t,R/1E16,c='k')
ax.text(0.05, 0.95, "$r \propto t^{1.54}$", fontsize=12, transform=ax.transAxes,
        ha='left', va='top')
ax.set_ylabel("Radius ($10^{16}$ cm)")

ax = axarr[2,0]
ax.plot(t,(R/t_sec)/3E10,c='k')
ax.text(0.05, 0.95, "$v \propto t^{0.54}$", fontsize=12, transform=ax.transAxes,
        ha='left', va='top')
ax.set_ylabel("Velocity ($c$)")
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.set_xlabel("Time (days)", fontsize=12)

ax = axarr[0,1]
ax.plot(R/1E16,B,c='k')
ax.text(0.95, 0.95, "$B \propto r^{-2.0}$", fontsize=12, transform=ax.transAxes,
        ha='right', va='top')
ax.set_ylabel("Magnetic Field ($G$)")

ax = axarr[1,1]
ax.plot(R/1E16,ne,c='k')
ax.set_ylabel("Electron Density (cm$^{-3}$)")
ax.text(0.95, 0.95, "$n_e \propto r^{-2.63}$", fontsize=12, transform=ax.transAxes,
        ha='right', va='top')

ax = axarr[2,1]
ax.plot(R/1E16, Mdot, c='k')
ax.set_xlabel("Radius ($10^{16}$) cm")
ax.set_ylabel("Mass-Loss Rate ($M_\odot$/yr)")
ax.text(0.95, 0.95, "$\dot{M} \propto r^{-0.63}$", fontsize=12, transform=ax.transAxes,
        ha='right', va='top')

for ax in axarr.flatten():
    ax.set_yscale('log')
    ax.set_xscale('log')

plt.tight_layout()
plt.show()
