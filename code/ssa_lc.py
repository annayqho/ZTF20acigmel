""" SSA light curves """
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from get_radio import *

alpha_r = 1
p = 3
s = 1.4

# CONSTANTS
c = 3E10
m_p = 1.67E-24
c6 = 8.16E-41
c5 = 6.29E-24
c1 = 6.27E18
El = 8.17E-7
f = 0.5
q_e = 4.8E-10
m_e = 9.1E-28
sigmat = 6.6524E-25
eps_e = 1/3
eps_B = 1/3


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


def Fnu(nu, t):
    """ calculate the expected light curve at a given frequency """

    # First: calculate the cooling frequency over time



    ft = np.zeros(len(t))

    choose = nu < nu_m(t)
    ft[choose] = Fpeak(t)[choose] \
            * (nu_m(t)[choose]/nu_a(t)[choose])**(5/2) \
            * (nu/nu_m(t)[choose])**2 

    choose = np.logical_and(nu > nu_m(t), nu < nu_a(t))
    ft[choose] = Fpeak(t)[choose] * (nu/nu_a(t)[choose])**(5/2)

    choose = nu > nu_a(t)
    ft[choose] = Fpeak(t)[choose] * (nu/nu_a(t)[choose])**(-(p-1)/2)
    return ft


def nu_m(t):
    return 2*(t/57)**((5-s/2)*alpha_r-5)


def nu_a(t):
    return 35*(t/57)**(((p*(5-s/2)-3*s)*alpha_r-5*p+2)/(p+4))


if __name__=="__main__":
    # Get the data
    islim, tel, freq, days, flux, eflux = get_data_all()

    # Get some fiducial parameters
    #run(57, 20, 0.5*1E-3, p, 0.2442)

    # Start with the 10 GHz light curve
    choose = np.logical_and(freq==12, islim==False)

    # Plot and take a look
    plt.errorbar(days[choose]/1.2442, flux[choose], eflux[choose], fmt='o', c='grey')
    
    t = np.linspace(5,100)

    nu = 12
    flux = Fnu(nu, t)

    plt.show()
