""" Fit a broken power law to the SED 

and determine the evolution of that SED with time
"""

from scipy.optimize import curve_fit
from get_radio import *

z = 0.2442
z = 0 

t0 = 72 / (1+z) # reference time


def func(x_in, Fa0, nua0, alpha1, alpha2, beta1, beta2, s):
    nu,t = x_in

    pref = Fa0 * (t/t0)**(alpha1)
    bot = nua0 * (t/t0)**(alpha2)
    Fnu = pref * ((nu/bot)**(-s*beta1) + (nu/bot)**(-s*beta2))**(-1/s)

    return Fnu


def func_no_s(x_in, Fa0, nua0, alpha1, alpha2, beta1, beta2):
    s = 100
    nu,t = x_in

    pref = Fa0 * (t/t0)**(alpha1)
    bot = nua0 * (t/t0)**(alpha2)
    Fnu = pref * ((nu/bot)**(-s*beta1) + (nu/bot)**(-s*beta2))**(-1/s)

    return Fnu


def func_no_spindex(x_in, Fa0, nua0, alpha1, alpha2, s):
    beta1 = 5/2
    beta2 = -2
    nu,t = x_in

    pref = Fa0 * (t/t0)**(alpha1)
    bot = nua0 * (t/t0)**(alpha2)
    Fnu = pref * ((nu/bot)**(-s*beta1) + (nu/bot)**(-s*beta2))**(-1/s)

    return Fnu


def func_spec(nu, Fa0, nua0, beta1, beta2, s):
    Fnu = Fa0 * ((nu/nua0)**(-s*beta1) + (nu/nua0)**(-s*beta2))**(-1/s)
    return Fnu


islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()
#keep = tel=='VLA'
#islim = islim[keep]
#tel = tel[keep]
#freq_obs = freq_obs[keep]
#days_obs = days_obs[keep]
#flux_obs = flux_obs[keep]
#eflux_obs = eflux_obs[keep]

# Put everything into the rest-frame before modeling
z = 0
freq = freq_obs[islim==False] * (1+z)
flux = flux_obs[islim==False] / (1+z)
eflux = eflux_obs[islim==False] / (1+z)
days = days_obs[islim==False] / (1+z)

# Times of observation
bins_obs = np.array([13, 18, 24, 28.3, 31.8, 38, 46, 51.9, 71, 95, 132])
bins = bins_obs / (1+z)

use_ind = np.arange(len(bins))[4:]
#use_ind = np.arange(len(bins))[0:8]

# Choose two bins for now
t = [] # time 
x = [] # frequency
y = [] # flux
ey = []
for ii in use_ind:
    bin = bins[ii]
    choose = np.abs(days-bin) < bin/10
    [t.append(bin) for i in np.arange(sum(choose))]
    [x.append(val) for val in freq[choose]]
    [y.append(val) for val in flux[choose]]
    [ey.append(val) for val in eflux[choose]]
t = np.array(t)
x = np.array(x)
y = np.array(y)
ey = np.array(ey)

ignore = np.logical_and(t>100, x<10)
t = t[~ignore]
x = x[~ignore]
y = y[~ignore]
ey = ey[~ignore]

print(t,x)

# Plot the data
#fig,ax = plt.subplots(1,1,figsize=(6,4))
fig,ax = plt.subplots(1,1,figsize=(4,3))

#col = ['#003f5c', '#bc5090', '#ffa600']
col = ['k', '#003f5c', '#2f4b7c', '#665191', '#a05195', '#d45087', '#f95d6a', '#ff7c43', '#ffa600']
j = 0
for ii in use_ind:
    choose = t==bins[ii]
    plt.errorbar(x[choose], y[choose], yerr=ey[choose], fmt='o', c=col[j])
    j += 1

# Now do a fit...
## Fa0, nua0, alpha1, alpha2, beta1, beta2, s
#p0 = np.array([0.55, 30, -2, -1, 1.5, -2, 2])
#p0 = np.array([0.55, 30, -2, -1, 1.5, -2])

# for mm only
p0 = np.array([0.55, 30, -2, -0.4, 1])

popt, pcov = curve_fit(
        func_no_spindex, (x, t), y, 
        p0=p0, sigma=ey, absolute_sigma=True, maxfev=10000)
        #bounds=([0.1,10,-4,-4,0.1,-4],[2,50,-0.1,-0.1,4,-0.1]))

for i in np.arange(len(p0)):
    print("%s +/- %s" %(popt[i], np.sqrt(pcov[i,i])))

nuplot = np.logspace(0.5,2.4)
j = 0
for ii in use_ind:
    tplot = np.array([bins[ii]]*len(nuplot))
    fplot = func_no_spindex((nuplot,tplot), *popt)
    plt.plot(nuplot, fplot, c=col[j], label=int(bins_obs[ii]))
    j += 1

# Calculate the chi squared
chisq = 0
for ii in use_ind:
    choose = t==bins[ii]
    f_model = func_no_spindex((x[choose],t[choose]), *popt)
    chisq += sum(((f_model - y[choose])/ey[choose])**2)
red_chisq = chisq / (len(t)-len(p0))
print("reduced X2: %s" %red_chisq)

plt.yscale('log')
plt.xscale('log')
plt.legend(loc='upper left', ncol=3)
plt.text(0.05, 0.7, r'$\chi^2=0.7$', transform=ax.transAxes, fontsize=12)
plt.ylim(1E-2, 1.5)
plt.xlabel("Freq. (GHz)", fontsize=14)
plt.ylabel("Flux density [mJy]", fontsize=14)
plt.tight_layout()
plt.show()
#plt.savefig("nup_Fp_fits_high_freq.png", dpi=200)
#plt.close()
