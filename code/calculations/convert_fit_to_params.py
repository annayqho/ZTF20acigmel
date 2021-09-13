""" Convert values from the fits to physical parameters """
import numpy as np
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/astronomy/papers_active/ZTF20acigmel/code")
from basics import get_z
z = get_z()
d_mpc = Planck15.luminosity_distance(z=z).value
d_cm = d_mpc * 3.086E24

# Constants
p = 3
eps_e = 1/3
eps_B = 1/3
eps = eps_e / eps_B
f = 0.5 # filling factor
c = 3E10
m_p = 1.67E-24
c6 = 8.16E-41
c5 = 6.29E-24
c1 = 6.27E18
El = 8.17E-7
q_e = 4.8E-10
m_e = 9.1E-28
sigmat = 6.6524E-25

# From the late-time data, we have
alpha_1 = -2.2
ealpha_1 = 0.1
alpha_2 = -0.88
ealpha_2 = 0.20
s = 1.0
e_s = 0.2
f_p = 0.68E-3
ef_p = 0.08E-3
nu_p = 22
enu_p = 1
t_d = 58

# For a constant shock, we have
alpha_1 = -2.1
ealpha_1 = 0.1
target_value = 1
alpha_2 = 9*alpha_1/19-target_value
ealpha_2 = ealpha_1 * 9/19
s = 0.78
e_s = 0.13
f_p = 0.83E-3
ef_p = 0.11E-3
nu_p = 23
enu_p = 1
t_d = 58
et_d = 1


print("First, estimate B and nu_c under two different assumptions.")

fast_cooling = False
slow_cooling = False

# Estimate the B and therefore nu_c under two different assumptions
print("Trying the fast-cooling regime, nup > nuc")
B = 0.14 * eps**(-4/13) * (f_p)**(-2/13) * (d_mpc)**(-4/13) * (nu_p/5)**(21/13) * (t_d/100)**(4/13)
nu_c = 18*np.pi*m_e*c*q_e / (sigmat**2 * B**3 * (t_d * 86400)**2)
if nu_c < nu_p: 
    fast_cooling=True
print("Trying the slow-cooling regime, nup < nuc")
B = 0.3*(eps_e)**(-4/19) * (f_p)**(-2/19) * (d_mpc)**(-4/19) * (nu_p/5)
nu_c = 18*np.pi*m_e*c*q_e / (sigmat**2 * B**3 * (t_d * 86400)**2)
if nu_c > nu_p: 
    slow_cooling=True

if np.logical_and(slow_cooling, fast_cooling):
    print("Something is wrong...how can it be both?")

elif np.logical_and(slow_cooling==True, fast_cooling==False):
    print("We are in the slow-cooling regime.")
    R = 5.1E15 * (eps_e)**(-1/19) * (f_p)**(9/19) * (d_mpc)**(18/19) * (nu_p/5)**(-1)
    eRdFa = 5.1E15 * (eps_e)**(-1/19) * (d_mpc)**(18/19) * (nu_p/5)**(-1) * (9/19) * f_p**(-10/19)
    eRdnua = 5.1E15 * (eps_e)**(-1/19) * (f_p)**(9/19) * (d_mpc)**(18/19) * (nu_p/5)**(-2) * (1/5)
    eR = np.sqrt((eRdFa*ef_p)**2 + (eRdnua*enu_p)**2)
    B = 0.3*(eps_e)**(-4/19) * (f_p)**(-2/19) * (d_mpc)**(-4/19) * (nu_p/5)
    eBdFa = 0.3*(eps_e)**(-4/19) * (-2/19) * (f_p)**(-21/19) * (d_mpc)**(-4/19) * (nu_p/5)
    eBdnua = 0.3*(eps_e)**(-4/19) * (f_p)**(-2/19) * (d_mpc)**(-4/19) * (1/5)
    eB = np.sqrt((eBdFa*ef_p)**2 + (eBdnua*enu_p)**2)

elif np.logical_and(slow_cooling==False, fast_cooling==True):
    print("We are in the fast-cooling regime.")
    R = 4.2E15 * eps**(-1/13) * (f_p)**(6/13) * (d_mpc)**(12/13) * (nu_p/5)**(-11/13) * (t/100)**(1/13)
    eRdFa = 4.2E15 * eps**(-1/13) * (d_mpc)**(12/13) * (nu_p/5)**(-11/13) * (t/100)**(1/13) * (6/13) * f_p**(-7/13)
    eRdnua = 4.2E15 * eps**(-1/13) * (f_p)**(6/13) * (d_mpc)**(12/13) * (t/100)**(1/13) * (11/13) * (nu_p/5)**(-24/13) * (1/5)
    eR = np.sqrt((eRdFa*eFa)**2 + (eRdnua*enu_p)**2)
    B = 0.14 * eps**(-4/13) * (f_p)**(-2/13) * (d_mpc)**(-4/13) * (nu_p/5)**(21/13) * (t/100)**(4/13)
    eBdFa = 0.14 * eps**(-4/13) * (-2/13) * (f_p)**(-15/13) * (d_mpc)**(-4/13) * (nu_p/5)**(21/13) * (t/100)**(4/13)
    eBdnua = 0.14 * eps**(-4/13) * (f_p)**(-2/13) * (d_mpc)**(-4/13) * (21/13) * (nu_p/5)**(8/13) * (1/5) * (t/100)**(4/13)
    eB = np.sqrt((eBdFa*eFa)**2 + (eBdnua*enu_p)**2)

else:
    print("Not sure what's going on here.")

# Now we have R and B in-hand
print("R: %s +/- %s x 10^16" %(np.round(R/1E16,2), np.round(eR/1E16,2)))
print("B: %s +/- %s" %(np.round(B,2), np.round(eB,2)))

# And the evolution is
print("Now estimating the temporal evolution of R and B")
R_exp = alpha_1*(9/19) + (alpha_2)*(-1)
B_exp = alpha_1*(-2/19) + (alpha_2)*(1)
n_exp = B_exp*2
print("R ~ t^(%s), B ~t^(%s), n~t^(%s)" %(np.round(R_exp,2), np.round(B_exp,2), np.round(n_exp,2)))

# Next, other things
print("Now we're going to measure other quantities")

v = R * (t_d*86400)**(-1)
dvdR = (t_d*86400)**(-1)
dvdt = R * (t_d*86400)**(-2) * (86400)**(-1)
ev = np.sqrt((dvdt*et_d)**2 + (dvdR*eR)**2)
print("shock speed %s +/- %s c" %(np.round(v/3E10,3), np.round(ev/3E10,3)))

#V = (4/3) * f * np.pi * R**3 # volume

P = B**2/(8*np.pi*eps_B)
eP = (1/(8*np.pi*eps_B)) * 2*B**1 * eB

n_e = (4*P/3) * (1/m_p) * v**(-2)
dnedP = (1/v**2) * (1/m_p) * (4/3)
dnedv = (4*P/3) * (1/m_p) * (-2)*v**(-3)
ene = np.sqrt((dnedP*eP)**2 + (dnedv*ev)**2)
print("ambient density %s +/- %s x 10^3 cm-3" %(
    np.round(n_e/1E3,2), np.round(ene/1E3,2)))

U = (1/eps_B) * (4*np.pi/3) * f * (R**3) * (B**2 / (8*np.pi))
dUdR = (1/eps_B) * (4*np.pi/3) * f * (B**2 / (8*np.pi)) * 3 * R**2
dUdB = (1/eps_B) * (4*np.pi/3) * f * (R**3) * (2*B / (8*np.pi))
eU = np.sqrt((dUdR*eR)**2 + (dUdB*eB)**2)
print("energy %s +/- %s x 10^48 cm-3" %(
    np.round(U/1E48,2), np.round(eU/1E48,2)))

v_w = 1000*3E5
Mdot = (n_e * 4 * np.pi * m_p * R**2 * v_w / 2E33)*(365*86400)
print("mass loss rate for 1000 km/s", Mdot)

v_w = 10*3E5
Mdot = (n_e * 4 * np.pi * m_p * R**2 * v_w / 2E33)*(365*86400)
print("mass loss rate for 10 km/s", Mdot)

gamma_m = 1 + eps_e * (m_p / m_e) * (v**2/c**2)
gamma_m = 1 + (1/2)*((p-2)/(p-1))*eps_e * (m_p / m_e) * (v**2/c**2)
print("gamma m", gamma_m)

nu_g = (q_e * B) / (2*np.pi*m_e*c)
nu_m = gamma_m**2 * nu_g
print("nu m", nu_m/1E9)

nu_c = 18*np.pi*m_e*c*q_e / (sigmat**2 * B**3 * (t_d * 86400)**2)
enu_c = 18*np.pi*m_e*c*q_e*B**(-4)*eB / (sigmat**2 * (t_d * 86400)**2)

print("nu c: %s +/- %s GHz" %(np.round(nu_c/1E9,2), np.round(enu_c/1E9,2)))
