""" Fit a relativistic Maxwellian """

import matplotlib.pyplot as plt
from get_radio import *
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15


def fitfunc(nu_ghz, f_m, tau_m, nu_m):
    """ Treat all frequencies in GHz """
    # Eq 29 from Mahadevan paper
    xM = 2*nu_ghz/(3*nu_m)
    Inu = 2.5651*(1+1.92/xM**(1/3)+0.9977/xM**(2/3))*np.exp(-1.8899*xM**(1/3))

    # Ben's equation
    exponent = -tau_m * (nu_ghz/nu_m)**(-1) * Inu
    fnu = f_m * xM**2 * (1-np.exp(exponent))

    return fnu


def fitfunc_physical(nu_ghz, B, ne):
    """ This time with physical parameters

    Te_rel = kT/mc^2
    beta = v/c
    R = in cm
    B = in Gauss

    Treat all frequencies in GHz """
    beta = 0.177

    # Scale factor in mJy
    v_cgs = beta*3E10
    L = 4.59E22 * B**2 * (beta/0.1)**12 * (td/50)**2
    dcm = Planck15.luminosity_distance(z=0.2442).cgs.value
    fmjy = (L / (4*np.pi*dcm**2)) / 1E-23 / 1E-3

    # Synchrotron frequency in GHz
    nu_m = (0.033) * (beta/0.1)**4 * B

    # taum
    tau_m = 1.18E6 * ne * B**(-1) * (beta/0.1)**(-9) * (td/50)

    # Eq 29 from Mahadevan paper
    xM = nu_ghz/nu_m
    Inu = 2.5651*(1+1.92/xM**(1/3)+0.9977/xM**(2/3))*np.exp(-1.8899*xM**(1/3))

    # Ben's equation
    exponent = -tau_m * (nu_ghz/nu_m)**(-1) * Inu
    fnu = fmjy * xM**2 * (1-np.exp(exponent))

    return fnu


def fitfunc_powlaw(x, A, beta):
    return A*x**(beta)


def fitfunc_exponential(nu, const, nuM):
    xM = nu/nuM
    return const*np.exp(-1.8899*xM**(1/3))


def at2020xnd(ax):
    islim, tel, freq, days, flux, eflux = get_data_all()
    z = 0.2442

    # Choose two epochs of observations to fit together
    choose = np.logical_and.reduce((days>41, days<52, islim==False))
    # Get values in rest-frame
    x = freq[choose] * (1+z)
    y = flux[choose] / (1+z)
    ey = eflux[choose] / (1+z)
    # Sort
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    ey = ey[order]
    td = np.average(days[choose]) / (1+z)

    # I think you can divide the x-value by (51.9/46)**1
    # to account for the differing dates

    # Plot the data
    col = '#a05195'
    marker = '*'
    msize = 14
    ax.errorbar(x, y, ey, 
            fmt='%s-' %marker, c=col, label=None, ms=msize)

    # Fit for a Maxwellian w/o physical parameters
    p0 = [0.0003, 5.5E4, 0.7]

    # Fit for the Maxwellian in terms of physical quantities
    popt, pcov = curve_fit(fitfunc, x, y, p0=p0, #maxfev=1000000,
            sigma=ey, absolute_sigma=True)
    xfit = np.linspace(1,300)
    yfit = fitfunc(xfit, *popt)
    print("Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))

    yfit = fitfunc(xfit, *p0)
    ax.plot(
            xfit,yfit, c='k', ls='--', zorder=5, 
            label=r'Maxwellian ($\nu_m\approx4$)')

    ax.set_xticks([10,30,50,100,200])
    ax.set_xticklabels([10,30,50,100,200])
    ax.set_yticks([0.1, 0.2, 0.3, 0.4, 0.6, 1])
    ax.set_yticklabels([0.1, 0.2, 0.3, 0.4, 0.6, 1])
    plt.minorticks_off()

    ax.set_xlim(9, 300)
    ax.set_ylim(7E-2, 1)

    #plt.savefig("camel_sed_maxwellian.png", dpi=300)
    #plt.close()


def at2020xnd_high_freq(ax):
    """ Fit only to the NOEMA data, and include the power law """

    islim, tel, freq, days, flux, eflux = get_data_all()
    z = 0.2442

    # Choose two epochs of observations to fit together
    choose = np.logical_and.reduce((days>41, days<52, islim==False))
    # Get values in rest-frame
    x = freq[choose] * (1+z)
    y = flux[choose] / (1+z)
    ey = eflux[choose] / (1+z)
    # Sort
    order = np.argsort(x)
    x = x[order][4:]
    y = y[order][4:]
    ey = ey[order][4:]
    td = np.average(days[choose]) / (1+z)

    # Plot the data
    col = '#a05195'
    marker = '*'
    msize = 14
    ax.errorbar(x, y, ey, 
            fmt='%s-' %marker, c=col, label=None, ms=msize)

    # Fit for a Maxwellian w/o physical parameters
    # This function is only the exponential, so a constant and nuM
    p0 = [0.0003, 0.7]
    popt, pcov = curve_fit(
            fitfunc_exponential, x, y, sigma=ey, absolute_sigma=True, p0=p0,
            bounds=((0, 0), (np.inf, np.inf)))

    #nsim = 300
    #ysamples = np.zeros((nsim, len(x)))
    #consts = np.zeros(nsim)
    #xms = np.zeros(nsim)
    #for jj,val in enumerate(y):
    #    ysamples[:,jj] = np.random.normal(loc=val, scale=ey[jj], size=nsim)
    #for jj in np.arange(nsim):
    #    popt, pcov = curve_fit(
    #            fitfunc_exponential, x, ysamples[jj], p0=p0,
    #            bounds=((0, 0), (np.inf, np.inf)))
    #    xmod = np.linspace(1,300)
    #    ymod = fitfunc_exponential(xmod, *popt)
    #    ax.plot(xmod, ymod, alpha=0.1, c='grey')
    #    consts[jj] = popt[0]
    #    xms[jj] = popt[1]

    xfit = np.linspace(1,300)
    yfit = fitfunc_exponential(xfit, *popt)
    print("Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
    ax.plot(xfit,yfit, c='k', ls='--', zorder=5, 
            label=r'Maxwellian ($\nu_m\approx4\,\mathrm{GHz}$)')

    # Fit for a power law
    popt, pcov = curve_fit(fitfunc_powlaw, x, y, maxfev=1000000)
    xfit = np.linspace(1,300)
    yfit = fitfunc_powlaw(xfit, *popt)
    print("Power law fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))

    # Plot a power law with p=3
    yfit = y[0]*(xfit/x[0])**(-1.5)
    ax.plot(xfit,yfit, c='k', ls='-', zorder=5, label='Power law ($p=3$)')



def at2018cow():
    fig,ax = plt.subplots(1,1,figsize=(5.5,4))

    x = np.array([9.0, 34.0, 243.3, 259.3, 341.5, 357.5])
    y = np.array([0.27, 5.6, 36.6, 31.21, 19.49, 17.42])
    ey = np.array([0.06, 0.16, 0.81, 0.92, 1.47, 2.8])
    td = 10.5

    # Plot the data
    col = '#a05195'
    marker = 'o'
    msize = 5
    ax.errorbar(x, y, ey, 
            fmt='%s-' %marker, c=col, label=None, ms=msize)

    # Fit for a Maxwellian w/o physical parameters
    p0 = [5E-2, 3E4, 1]
    popt, pcov = curve_fit(fitfunc, x, y, sigma=ey, absolute_sigma=True, p0=p0, maxfev=1000000)
    xfit = np.linspace(1,400)
    yfit = fitfunc(xfit, *popt)
    print("Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
    ax.plot(xfit,yfit, c='k', ls='--', zorder=5, label='Maxwellian')

    ax.set_xscale('log')
    ax.set_yscale('log')

    #ax.set_xticks([10,30,50,100,200])
    #ax.set_xticklabels([10,30,50,100,200])
    #ax.set_yticks([0.2, 1, 3, 5, 10, 50])
    #ax.set_yticklabels([0.2, 1, 3, 5, 10, 50])
    plt.minorticks_off()

    ax.set_xlabel("Rest Frequency [GHz]", fontsize=14)
    ax.set_ylabel("Rest Flux Density (mJy)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    ax.set_xlim(7, 400)
    ax.set_ylim(0.18, 50)

    plt.tight_layout()
    plt.show()
    #plt.savefig("cow_sed_maxwellian.png", dpi=300)


def ultralong():
    fig,ax = plt.subplots(1,1,figsize=(5.5,4))

    z = 0.347
    x = np.array([4.8, 7.4, 9.5, 13.5, 16.0, 22.0]) * (1+z)
    y = np.array([237, 298, 293, 146, 89, 104]) / (1+z)
    ey = np.array([17, 12, 14, 20, 22, 13]) / (1+z)
    td = 2

    # Plot the data
    col = '#a05195'
    marker = 'o'
    msize = 5
    ax.errorbar(x, y, ey, 
            fmt='%s-' %marker, c=col, label=None, ms=msize)

    # Fit for a Maxwellian w/o physical parameters
    p0 = [5E-1, 1000, 1E-1]

    nsim = 100
    ysamples = np.zeros((nsim, len(x)))
    nums = np.zeros(nsim)
    taums = np.zeros(nsim)
    fms = np.zeros(nsim)
    for jj,val in enumerate(y):
        ysamples[:,jj] = np.random.normal(loc=val, scale=ey[jj], size=nsim)
    for jj in np.arange(nsim):
        popt, pcov = curve_fit(fitfunc, x, y, p0=p0, maxfev=1000000)
        fms[jj] = popt[0]
        taums[jj] = popt[1]
        nums[jj] = popt[2]
    xfit = np.linspace(1,30)
    yfit = fitfunc(xfit, *popt)
    print("Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
    ax.plot(xfit,yfit, c='k', ls='--', zorder=5, label='Maxwellian')

    plt.tight_layout()
    plt.show()
    #plt.savefig("cow_sed_maxwellian.png", dpi=300)


if __name__=="__main__":
    # Initialize figure
    fig,axarr = plt.subplots(2,1,figsize=(4,6), sharey=True)

    # Top panel
    at2020xnd_high_freq(axarr[0])

    # Bottom panel
    at2020xnd(axarr[1])

    # Formatting
    for ax in axarr:
        ax.set_ylabel("Rest Flux Density (mJy)", fontsize=14)
        ax.tick_params(axis='both', labelsize=12)
        ax.set_ylim(7E-2, 1)
        ax.set_yticks([0.1,0.2,0.3,0.6,1])
        ax.set_yticklabels([0.1,0.2,0.3,0.6,1])
    axarr[1].set_xlabel("Rest Frequency [GHz]", fontsize=14)
    ax = axarr[0]
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim(90, 300)
    ax.set_xticks([100,150,200,300])
    ax.set_xticklabels([100,150,200,300])
    ax = axarr[1]
    ax.scatter(10*1.2442, 0.124/1.2442)
    ax.set_xlim(7, 400)
    ax.set_xticks([10, 30, 100, 300])
    ax.set_xticklabels([10, 30, 100, 300])

    plt.minorticks_off()
    
    # Display or save
    plt.tight_layout()
    plt.show()
