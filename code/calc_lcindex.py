""" Calculate the index of the decay of the peak flux with frequency """
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import warnings


def func(x,a,beta):
    return a*(x**beta)


def fit():
    xvals = np.array([94,79,36,26.5, 8])
    yvals = np.array([1.03, 1.09, 0.68, 0.50, 0.18])
    eyvals = np.array([0.04, 0.05, 0.17, 0.01, 0.02])
    plt.errorbar(xvals, yvals, eyvals, fmt='o', c='k')

    xvals = xvals[2:]
    yvals = yvals[2:]
    eyvals = eyvals[2:]
    plt.errorbar(xvals, yvals, eyvals, fmt='o', c='k')

    nsim = 600

    ysamples = np.zeros((nsim, len(xvals)))
    for jj,val in enumerate(yvals):
        ysamples[:,jj] = np.random.normal(
                loc=val,scale=eyvals[jj],size=nsim)

    betas = np.zeros(nsim)
    alphas = np.zeros(nsim)

    for jj in np.arange(nsim):
        plt.scatter(xvals,ysamples[jj],marker='.',c='k',zorder=5)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            popt, pcov = curve_fit(
                    func, xvals, ysamples[jj], p0=[1,1])
                    #bounds=([1000,0],[50000, np.inf]),maxfev=100000)
            betas[jj] = popt[1]
            alphas[jj] = popt[0]
            xplot = np.linspace(1,100)
            yplot = func(xplot, popt[0], popt[1])
            plt.plot(xplot,yplot,lw=0.5,alpha=0.5,c='lightgrey')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(8,100)

def plot():
    fig = plt.figure(figsize=(4,3))
    xvals = np.array([94,79,36,26.5, 8])
    yvals = np.array([1.03, 1.09, 0.68, 0.50, 0.18])
    eyvals = np.array([0.04, 0.05, 0.17, 0.01, 0.02])
    plt.errorbar(xvals, yvals, eyvals, fmt='o', c='k')
    xplot = np.linspace(1,100)
    yplot = func(xplot, 0.07, 0.62)
    plt.plot(
            xplot, yplot, lw=0.5, c='k', 
            label=r"$F_{\mathrm{pk},\nu} \propto \nu^{0.6}$")
    xplot = np.linspace(1,100)
    yplot = func(xplot, 0.032, 0.85)
    plt.plot(xplot, yplot, lw=0.5, c='k', ls='--',
            label=r"$F_{\mathrm{pk},\nu} \propto \nu^{0.8}$")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(6,110)
    plt.ylim(0.1,1.5)
    plt.legend()
    plt.xlabel("Observing Frequency [GHz]", fontsize=12)
    plt.ylabel("Peak Flux Density", fontsize=12)
    plt.tick_params(axis='both', labelsize=12)
    plt.tight_layout()
    plt.savefig("fpeak_nu.png", dpi=200)
    plt.close()
    #plt.show()


if __name__=="__main__":
    plot()
