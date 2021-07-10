""" SSA light curves """
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from get_radio import *


# Explosion Constants
alpha_r = 1.0
p = 3
eps_e = 1/3
eps_B = 1/3
eps = eps_e / eps_B
f = 0.5 # filling factor
z = 0.2442
d_mpc = Planck15.luminosity_distance(z=z).value
d_cm = Planck15.luminosity_distance(z=z).cgs.value

# Fiducial Values
t0 = 72/(1+z)

# Frequencies to fit for
#use_nu = np.array([117, 98])/(1+z)
#use_nu = np.array([45,33])/(1+z)
use_nu = np.array([117,98,45,33,22,15,10])/(1+z)
#use_nu = np.array([33,22,15,10])/(1+z)
#use_nu = [15]

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


def Fnu(nu, t, k, s, opt_thick_index, Fa0, nua0):
    """ calculate the expected light curve at a given frequency 

    all quantities are in the rest-frame
    """
    t_sec = t*86400 # original t is in days, need it in sec for cgs things

    t_model = np.logspace(1,3,1000)
    t_model_sec = t_model*86400

    # Step 1: Calculate the magnetic field strength
    # this equation assumes nu_a < nu_c, valid at the late times (72d) considered here
    B0 = 0.30*(eps)**(-4/19)*(Fa0)**(-2/19)*(d_mpc)**(-4/19)*(nua0/5)
    B = B0*(t_model/t0)**(-alpha_r*k/2) # model time, not observed times

    # Step 2: Calculate the cooling frequency
    nu_c = 18*np.pi*m_e*c*q_e / (sigmat**2 * B**3 * t_model_sec**2) # model time

    # Step 3: Assert nu_c < nu_a (as at early times), and check
    exponent = alpha_r-1-alpha_r*(k*(p+3))/(2*(p+5))
    nua1 = nua0 * (t_model/t0)**exponent # model time
    exponent = 4*alpha_r-2-alpha_r*(k*(2*p+5)/(2*(p+5)))
    Fa1 = Fa0 * (t_model/t0)**exponent # model time
    is_right_1 = nu_c/1E9 < nua1

    # Step 4: Using nu_c = nu_a from above as the intersection point,
    # use the expressions for nu_c > nu_a afterwards
    top_exp = 2*(p+6)-2*alpha_r*(p+8)+alpha_r*k*(p+6)
    bot_exp = 2*(p+4)
    nu_exponent = -top_exp/bot_exp
    f_exponent = -((2*p+13)*(2-alpha_r*(4-k))) / (2*(p+4))
    if sum(is_right_1)>0:
        nua2 = nua1[is_right_1][-1] * (t_model/t_model[is_right_1][-1])**exponent
        Fa2 = Fa1[is_right_1][-1] * (t_model/t_model[is_right_1][-1])**exponent
    else:
        nua2 = nua0 * (t_model/t0)**exponent
        Fa2 = Fa0 * (t_model/t0)**exponent
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

    F_res = np.interp(t,t_model,F)
    is_right_res = np.interp(t,t_model,is_right_1)
    nua_res = np.interp(t,t_model,nua)
    Fa_res = np.interp(t,t_model,Fa)

    return F_res, is_right_res, nua_res, Fa_res


def nu_m(t):
    return 2*(t/57)**((5-s/2)*alpha_r-5)


def nu_a(t):
    return 35*(t/57)**(((p*(5-s/2)-3*s)*alpha_r-5*p+2)/(p+4))


def fit_func(x_in, k, s, opt_thick_index, Fa0, nua0):
    """ The fitting function """
    t,nplt = x_in
    out = []
    for nu in use_nu:
        choose = nplt==nu
        order = np.argsort(t[choose])
        Fnut = Fnu(nu, t[choose][order], k, s, opt_thick_index, Fa0, nua0)[0]
        [out.append(val) for val in Fnut]
    out = np.array(out)

    return out


if __name__=="__main__":
    # Fit Parameters
    k = 2.7
    s = 1 # spectral smoothness
    opt_thick_index = 1.5
    Fa0 = 0.4 # mJy
    nua0 = 18 # GHz
    p0 = [k, s, opt_thick_index, Fa0, nua0]

    islim, tel, freq_obs, days, flux, eflux = get_data_all()
    freq = freq_obs / (1+z) # also need this in rest-frame

    t = []
    y = []
    ey = []
    nplt = []

    for nu in use_nu:
        if np.logical_or(nu<26, nu>28):
            choose = np.logical_and(freq==nu, islim==False)
        else:
            freq_crit = np.logical_and(freq>26, freq<28)
            choose = np.logical_and(freq_crit, islim==False)
        [t.append(val) for val in list(days[choose])]
        [y.append(val) for val in list(flux[choose])]
        [ey.append(val) for val in list(eflux[choose])]
        [nplt.append(val) for val in [nu]*sum(choose)]
    t_obs = np.array(t)
    y_obs = np.array(y)
    ey_obs = np.array(ey)
    nplt = np.array(nplt)

    # Need to convert t to rest-frame. frequencies are already rest-frame
    t = t_obs/(1+z)
    y = y_obs/(1+z)
    ey = ey_obs/(1+z) 

    # Plot the data
    cols = ['k', 'darkgrey', '#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600', 'lightgrey']
    #cols = ['darkgrey', 'lightgrey']

    for ii,nu in enumerate(use_nu):
        choose = nplt==nu
        plt.errorbar(
                t[choose], y[choose], ey[choose], fmt='o', 
                c=cols[ii], label=r"$%s$ GHz" %int(nu*(1+z)))

    # Fit for the parameters
    #popt, pcov = curve_fit(
    #        fit_func, (t,nplt), y, p0=p0, sigma=ey, absolute_sigma=True,
    #        bounds=([1,0.5,0.5,0.2,10],[4,5,3,2,50]))
    popt, pcov = curve_fit(
            fit_func, (t,nplt), y, p0=p0, sigma=ey, absolute_sigma=True,
            bounds=([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf],
                [np.inf,np.inf,np.inf,np.inf,np.inf]))
    print(popt)
    for ii in np.arange(len(popt)):
        print(np.sqrt(pcov[ii,ii]))

    # Plot the best-fit curve
    plot_t = np.logspace(-1,3,1000)
    #popt = np.array([2.7, 1, 1.6, 1, 18])
    for ii,nu in enumerate(use_nu):
        out = Fnu(nu, plot_t, *popt)[0]
        plt.plot(plot_t[out>0], out[out>0], c=cols[ii])

    out = Fnu(nu, plot_t, *popt)

    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(0.01,2)
    plt.xlim(8,200)
    plt.legend(loc='upper right', fontsize=9, ncol=4)
    plt.xlabel("Days", fontsize=16)
    plt.ylabel("Flux Density (mJy)", fontsize=16)
    plt.tick_params(axis='both', labelsize=14)
    plt.tight_layout()
    plt.show()
    #plt.savefig("model_fits.png", dpi=200)
    #plt.close()
