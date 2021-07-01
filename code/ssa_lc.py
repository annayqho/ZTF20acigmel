""" SSA light curves """
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from get_radio import *

use_nu = [94, 79, 26, 12, 8]

# Explosion Constants
p = 3
eps_e = 0.1
eps_B = 0.1
eps = eps_e / eps_B
f = 0.5 # filling factor
z = 0.2442
d_mpc = Planck15.luminosity_distance(z=z).value
d_cm = Planck15.luminosity_distance(z=z).cgs.value


# Numerical Constants
c = 3E10
m_p = 1.67E-24
c6 = 8.16E-41
c5 = 6.29E-24
c1 = 6.27E18
El = 8.17E-7
q_e = 4.8E-10
m_e = 9.1E-28
sigmat = 6.6524E-25


def get_R(Fp, nup, d_mpc):
    """ This is Equation 13 in C98
    Fp in Jy
    nup in GHz
    Assumes that the power-law index gamma = 3

    Returns Rp in cm
    """
    alpha = eps_e / eps_B
    Rp = 8.8E15 * alpha**(-1/19) * (f/0.5)**(-1/19) * (Fp)**(9/19) * \
         d_mpc**(18/19) * (nup/5)**(-1)
    return Rp


def get_B(Fp, nup, d_mpc):
    """ This is Equation 14 in C98
    Fp in Jy
    nup in GHz

    Returns Bp in Gauss
    """
    alpha = eps_e / eps_B
    Bp = 0.58 * alpha**(-4/19) * (f/0.5)**(-4/19) * (Fp)**(-2/19) * \
         d_mpc**(-4/19) * (nup/5)
    return Bp


def get_U(Fp, nup, d_mpc):
    """ This is Equation 12 in our paper
    Fp in Jy, nup in GHz, D is the ANGULAR DIAMETER DISTANCE in Mpc """
    prefactor = 1.92E46
    U = prefactor * (1/eps_B) * (eps_e/eps_B)**(-11/19) * (f/0.5)**(8/19) * \
            (Fp)**(23/19) * (d_mpc)**(46/19) * (nup/5)**(-1)
    return U


def run(dt, nupeak, fpeak, p, z):
    # To convert from flux density to luminosity,
    # you need to use the luminosity distance
    DL_cm = Planck15.luminosity_distance(z=z).cgs.value
    print("lpeak", fpeak*1E-23*4*np.pi*DL_cm**2)

    # In the Chevalier 98 equation, D is the angular diameter distance
    DA_mpc = Planck15.angular_diameter_distance(z=z).value
    R = get_R(fpeak, nupeak, DA_mpc)
    print("R", R/10**15)
    V = (4/3) * f * np.pi * R**3 # volume

    B = get_B(fpeak, nupeak, DA_mpc)
    print("B", B)

    # Convert time to rest-frame...or maybe not?
    dt_rest = dt/(1+z)
    #dt_rest = dt
    v = R/(86400*dt_rest)
    print("v/c", v/3E10)
    gammabeta = v/c
    print("gamma beta", gammabeta)
    beta = 1/(1+(3E10/v)**2)**0.5
    print("Beta", beta)
    gamma = 1/(1-beta**2)
    print("gamma", gamma)

    P = B**2/(8*np.pi*eps_B)
    rho = (4*P/3)/v**2
    print("rho", rho)
    n_p = rho/m_p
    n_e = n_p
    print("ne", n_e)
    UB = V * (B**2)/(8*np.pi)

    U = get_U(fpeak, nupeak, DA_mpc)
    print("E", U)
    tauff = 8.235E-2 * n_e**2 * (R/3.086E18) * (8000)**(-1.35)
    print("tau_ff", tauff)
    Te = 1.2E6
    nuff = (1/(68*(Te/8000)**(-1.35)))**(-1/2.1)
    print("nu_ff [GHz]", nuff)
    Lff = 1.43E-27 * n_e**2 * Te**(1/2) * (4/3) * np.pi * (6E16)**3
    print("L_ff", Lff)
    # Gamma_m
    gammam = eps_e * (m_p / m_e) * ((p-2)/(p-1)) * beta**2
    print("gamma_m", gammam)
    # Gyrofrequency
    gammag = q_e*B / (2 * np.pi * m_e * c)
    print("gyrofrequency [MHz]", gammag/1E6)
    # nu_m
    num = gammam**2 * gammag
    print("nu_m [GHz]", num/1E9)
    # Cooling frequency
    gammac = 6 * np.pi * m_e * c / (sigmat * B**2 * dt * 86400)
    print("gammac", gammac)
    nuc = gammac**2 * q_e * B / (2 * np.pi * m_e * c)
    print("nu_cc [GHz]", nuc/1E9)


def Fpeak(t):
    top = (p*(5-s/2)-3*s)*alpha_r - 5*p + 2
    bottom = (p+4) * ((5-s/2)*alpha_r - 5)
    exponent = (4-3*s/2)*alpha_r - 1 - ((p-1)/2)*(top/bottom)
    return 0.5*(t/57)**exponent


def Fnu(nu, t, alpha_r, k, s, opt_thick_index, t0, Fa0, nua0):
    """ calculate the expected light curve at a given frequency """
    t_sec = t*86400 # original t is in days, need it in sec for cgs things

    # Step 1: Calculate the magnetic field strength
    # this equation assumes nu_a > nu_c, probably valid at early times
    B0 = 0.14*(eps)**(-4/13)*(Fa0)**(-2/13)*(d_mpc)**(-4/13)*(nua0/5)**(21/13)*(t0/100)**(4/13)
    B = B0*(t/t0)**(-alpha_r*k/2)

    # Step 2: Calculate the cooling frequency
    nu_c = 18*np.pi*m_e*c*q_e / (sigmat**2 * B**3 * t_sec**2)

    # Step 3: Assert nu_c < nu_a (as at early times), and check
    exponent = 2*(alpha_r-1)/(2*p+5) - ((p+3)*alpha_r*k)/(2*(p+5))
    nua1 = nua0 * (t/t0)**exponent
    exponent = ((2*p+15)*alpha_r - (2*p+5)*alpha_r*k/2 - 5) / (p+5)
    Fa1 = Fa0 * (t/t0)**exponent
    is_right_1 = nu_c/1E9 < nua1

    # Step 4: Assert nu_c > nu_a (as at late times) and check
    exponent = alpha_r*((4-k*(p+6))/(2*(p+4)))
    nua2 = nua0 * (t/t0)**exponent
    exponent = alpha_r*(2*p+13)*(2-k)/(2*(p+4))
    Fa2 = Fa0 * (t/t0)**exponent
    is_right_2 = nu_c/1E9 > nua2

    # Step 5: Solve for Fa(t)
    Fa = np.zeros(len(nu_c))
    Fa[is_right_1] = Fa1[is_right_1]
    Fa[is_right_2] = Fa2[is_right_2]
    nua = np.zeros(len(nu_c))
    nua[is_right_1] = nua1[is_right_1]
    nua[is_right_2] = nua2[is_right_2]

    # Step 6: Solve for F(nu,t)
    F = np.zeros(len(nu_c))

    choose = is_right_1
    F[choose] = Fa[choose]*((nu/nua[choose])**(-s*(opt_thick_index))+\
            (nu/nua[choose])**(-s*(-p/2)))**(-1/s)

    choose = is_right_2
    F[choose] = Fa[choose]*((nu/nua[choose])**(-s*(opt_thick_index))+\
            (nu/nua[choose])**(-s*(-(p-1)/2)))**(-1/s)

    return F,is_right_1,nua,Fa


def nu_m(t):
    return 2*(t/57)**((5-s/2)*alpha_r-5)


def nu_a(t):
    return 35*(t/57)**(((p*(5-s/2)-3*s)*alpha_r-5*p+2)/(p+4))


def fit_func(x_in, alpha_r, k, s, opt_thick_index, t0, Fa0, nua0):
    """ The fitting function """
    alpha_r=1
    t,nplt = x_in
    out = []
    for nu in use_nu:
        choose = nplt==nu
        Fnut = Fnu(nu, t[choose], alpha_r, k, s, opt_thick_index, t0, Fa0, nua0)[0]
        [out.append(val) for val in Fnut]
    out = np.array(out)

    return out


if __name__=="__main__":
    # Fit Parameters
    alpha_r = 1
    k = 3.8
    s = 1 # spectral smoothness
    opt_thick_index = 1.2
    t0 = 25 # days
    Fa0 = 3 # mJy
    nua0 = 110 # GHz
    p0 = [alpha_r, k, s, opt_thick_index, t0, Fa0, nua0]

    islim, tel, freq, days, flux, eflux = get_data_all()
    t = []
    y = []
    ey = []
    nplt = []

    for nu in use_nu:
        if nu!=26:
            choose = np.logical_and(freq==nu, islim==False)
        else:
            freq_crit = np.logical_and(freq>25, freq<28)
            choose = np.logical_and(freq_crit, islim==False)
        [t.append(val) for val in list(days[choose]/1.2442)]
        [y.append(val) for val in list(flux[choose])]
        [ey.append(val) for val in list(eflux[choose])]
        [nplt.append(val) for val in [nu]*sum(choose)]
    t = np.array(t)
    y = np.array(y)
    ey = np.array(ey)
    nplt = np.array(nplt)

    # Plot the data
    cols = ['#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600']

    for ii,nu in enumerate(use_nu):
        choose = nplt==nu
        plt.errorbar(t[choose], y[choose], ey[choose], fmt='o', c=cols[ii], label=r"$%s$ GHz" %nu)

    # Fit for the parameters
    popt, pcov = curve_fit(fit_func, (t,nplt), y, p0=p0)
    print(popt)
    for ii in np.arange(len(popt)):
        print(np.sqrt(pcov[ii,ii]))

    # Plot the best-fit curve
    plot_t = np.logspace(1,3,1000)
    for ii,nu in enumerate(use_nu):
        out = Fnu(nu, plot_t, *popt)[0]
        plt.plot(plot_t[out>0], out[out>0], c=cols[ii])

    out = Fnu(nu, plot_t, *popt)

    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(0.01,2)
    plt.xlim(8,200)
    plt.legend(loc='upper right', fontsize=12)
    plt.xlabel("Days", fontsize=16)
    plt.ylabel("Flux Density (mJy)", fontsize=16)
    plt.tick_params(axis='both', labelsize=14)
    plt.tight_layout()
    plt.show()
