""" Plot the spectral index (or just spectrum) of the ultra-long GRB130925A """
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def fitfunc(x,a,beta):
    return a*(x**beta)

dt = 2 # days
z = 0.347
freq = np.array([4.8, 7.4, 9.5, 13.5, 16.0, 22.0])
f = np.array([237, 298, 293, 146, 89, 104])
ef = np.array([17, 12, 14, 20, 22, 13])
plt.errorbar(freq, f, ef, fmt='o', c='lightgrey')

#freq = freq * (1+z)
#f = f / (1+z)
#ef = ef / (1+z)
#plt.errorbar(freq, f, ef, fmt='o', c='k')

plt.errorbar(freq[1:], f[1:], ef[1:], fmt='o', c='k')
popt,pcov = curve_fit(fitfunc, freq[1:], f[1:], sigma=ef[1:], absolute_sigma=False)
print(popt[1], np.sqrt(pcov[1,1]))
xplot = np.linspace(5,23)
yplot = fitfunc(xplot,*popt)
plt.plot(xplot, yplot, c='k', label=r'$f_\nu \propto \nu^{-1}$')

plt.errorbar(freq[2:], f[2:], ef[2:], fmt='o', c='k', mec='red')
popt,pcov = curve_fit(fitfunc, freq[2:], f[2:], sigma=ef[2:], absolute_sigma=False)
print(popt[1], np.sqrt(pcov[1,1]))
xplot = np.linspace(5,23)
yplot = fitfunc(xplot,*popt)
plt.plot(xplot, yplot, c='red', ls='--', label=r'$f_\nu \propto \nu^{-1.5}$')

plt.legend(fontsize=12)
plt.xlabel("Frequency", fontsize=13)
plt.ylabel("Flux density", fontsize=13)
plt.ylim(31, 359)
plt.tight_layout()
plt.show()

#plt.savefig("grb130925a.png", dpi=200)
#plt.close()
