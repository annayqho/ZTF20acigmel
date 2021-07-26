""" Fit a relativistic Maxwellian """

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset,inset_axes
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
    xM = 2*nu/(3*nuM)
    Inu = (1+1.92/xM**(1/3)+0.9977/xM**(2/3))*np.exp(-1.8899*xM**(1/3))
    #return const*nu*np.exp(-1.8899*xM**(1/3))
    return const*nu*Inu


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
    popt, pcov = curve_fit(fitfunc, x[1:], y[1:], p0=p0, #maxfev=1000000,
            sigma=ey[1:], absolute_sigma=True, 
            bounds=((0.0001,1E4,0.1),(0.0005,7E4,1.2)))
    xfit = np.linspace(1,300)
    yfit = fitfunc(xfit, *popt)
    print("Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))

    #yfit = fitfunc(xfit, *p0)
    yfit = fitfunc(xfit, *popt)
    ax.plot(
            xfit,yfit, c='k', ls='--', zorder=5, 
            label=r'Maxwellian ($\nu_m\approx4$)')

    ax.set_xticks([50,100,150,200,250])
    ax.set_xticklabels([50,100,150,200,250])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8])
    ax.set_yticklabels([0.2, 0.4, 0.6, 0.8])
    plt.minorticks_off()
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 0.9)



def at2020xnd_low_freq_late(ax):
    islim, tel, freq, days, flux, eflux = get_data_all()
    z = 0.2442

    bins = [71,95,132]
    col = ['#003f5c', '#bc5090', '#ffa600']

    for i,b in enumerate(bins):
        choose = np.logical_and.reduce((days>b-3, days<b+3, islim==False))
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

        if b==132:
            x = x[1:]
            y = y[1:]
            ey = ey[1:]

        # Plot the data
        marker = 'o'
        msize = 5
        ax.errorbar(x, y, ey, 
                fmt='%s-' %marker, c=col[i], ms=msize, label=str(b))

        # Fit for a Maxwellian w/o physical parameters
        p0 = [0.0003, 5.5E4, 0.7]

        # Fit for the Maxwellian in terms of physical quantities
        popt, pcov = curve_fit(fitfunc, x, y, p0=p0, #maxfev=1000000,
                sigma=ey, absolute_sigma=True, 
                bounds=((0.00001,0.01E4,0.1),(0.0006,6E4,1.0)))
        xfit = np.linspace(1,100)
        yfit = fitfunc(xfit, *popt)
        print("Maxwellian fit:")
        for i,param in enumerate(popt):
            print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))

        yfit = fitfunc(xfit, *popt)
        ax.plot(
                xfit,yfit, c='k', ls='--', zorder=5)

        ax.set_xticks([10,20,30,40,50])
        ax.set_xticklabels([10,20,30,40,50])
        ax.set_yticks([0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([0.2, 0.4, 0.6, 0.8])
        plt.minorticks_off()
        ax.set_xlim(7, 60)
        ax.set_ylim(0, 0.45)
        ax.legend()


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
    x = x[order][1:]
    y = y[order][1:]
    ey = ey[order][1:]
    td = np.average(days[choose]) / (1+z)

    # Plot the data
    col = '#a05195'
    marker = '*'
    msize = 14
    ax.errorbar(x, y, ey, 
            fmt='%s-' %marker, c=col, label=None, ms=msize)

    # Fit for a Maxwellian w/o physical parameters
    # This function is only the exponential, so a constant and nuM
    #p0 = [0.0003, 0.7]
    p0 = [0.0003, 100]
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
            label=r'Maxwellian ($\nu_m\approx1\,\mathrm{GHz}$)')

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
    fig,ax = plt.subplots(1,1,figsize=(5,4))

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

    axins = inset_axes(ax, 1.5, 1.5, loc=4, bbox_to_anchor=(0.9,0.25),
            bbox_transform=ax.figure.transFigure)
    axins.errorbar(x[2:], y[2:], ey[2:],fmt='%s-' %marker, c=col, 
            label=None, ms=msize)
    xfit = np.linspace(230,380)
    yfit = fitfunc(xfit, *popt)
    axins.plot(xfit,yfit, c='k', ls='--', zorder=5, label='Maxwellian')
    mark_inset(ax, axins, loc1=2, loc2=4)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xticks([10,20,30,50,100,200,300])
    ax.set_xticklabels([10,20,30,50,100,200,300])
    ax.set_yticks([0.2,0.4, 1, 3, 10, 30, 50])
    ax.set_yticklabels([0.2,0.4, 1, 3, 10, 30, 50])
    plt.minorticks_off()

    ax.set_xlabel("Rest Frequency [GHz]", fontsize=14)
    ax.set_ylabel("Rest Flux Density (mJy)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    ax.set_xlim(7, 400)
    ax.set_ylim(0.18, 80)

    plt.tight_layout()
    plt.show()
    #plt.savefig("cow_sed_maxwellian.png", dpi=300)
    #plt.close()


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
    #p0 = [5E-1, 1000, 1E-1]

    # just the exponential
    p0 = [5E-1, 1]

    #nsim = 100
    #ysamples = np.zeros((nsim, len(x)))
    #nums = np.zeros(nsim)
    #taums = np.zeros(nsim)
    #fms = np.zeros(nsim)
    #for jj,val in enumerate(y):
    #    ysamples[:,jj] = np.random.normal(loc=val, scale=ey[jj], size=nsim)
    #for jj in np.arange(nsim):
    #    #popt, pcov = curve_fit(fitfunc, x, y, p0=p0, maxfev=1000000)
    #    popt, pcov = curve_fit(fitfunc, x, y, p0=p0, maxfev=1000000)
    #    fms[jj] = popt[0]
    #    taums[jj] = popt[1]
    #    nums[jj] = popt[2]

    popt, pcov = curve_fit(
            fitfunc_exponential, x[2:5], y[2:5], sigma=ey[2:5], absolute_sigma=True, 
            p0=p0, maxfev=100000)
    xfit = np.linspace(1,30)
    yfit = fitfunc_exponential(xfit, *popt)
    print("Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
    ax.plot(xfit,yfit, c='k', ls='--', zorder=5, label='Maxwellian')

    xplot = np.linspace(1,20)
    yplot = (y[1]-10)*(xplot/x[1])**(1/3)
    ax.plot(xplot, yplot,c='k',ls='-',label=r'$f_{\nu} \propto \nu^{1/3}$')

    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(45,250)
    ax.set_xlim(5,34)
    ax.set_yticks([50,70,100,150,200])
    ax.set_yticklabels([50,70,100,150,200])
    ax.set_xticks([5,7,10,15,20,30])
    ax.set_xticklabels([5,7,10,15,20,30])
    ax.set_ylabel("Rest Flux Density (mJy)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel("Rest Frequency [GHz]", fontsize=14)
    plt.minorticks_off()
    plt.tight_layout()
    #plt.show()
    plt.savefig("ultralong_maxwellian.png", dpi=300)
    plt.close()


def make_camel_plots():
    """ Generate plots of the Camel """

    # two-panel
    #fig,axarr = plt.subplots(2,1,figsize=(4,6), sharey=True)

    # one-panel
    fig,ax = plt.subplots(1,1,figsize=(4,3), sharey=True)

    # Top panel
    #ax = axarr[0]
    #at2020xnd_high_freq(ax)

    # Bottom panel
    #ax = axarr[1]
    at2020xnd_low_freq_late(ax)

    # Formatting
    #for ax in axarr:
    ax.set_ylabel("Rest Flux Density (mJy)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)
    #ax.set_ylim(7E-2, 1)
    #ax = axarr[1]
    ax.set_xlabel("Rest Frequency [GHz]", fontsize=14)
    #ax = axarr[0]
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #ax.set_yticks([0.1,0.2,0.3,0.6,1])
    #ax.set_yticklabels([0.1,0.2,0.3,0.6,1])
    #ax.legend(loc='upper right', fontsize=9)
    #ax.set_xlim(90, 300)
    #ax.set_xticks([100,150,200,300])
    #ax.set_xticklabels([100,150,200,300])
    #ax = axarr[1]
    #ax.scatter(10*1.2442, 0.124/1.2442)
    #ax.set_xlim(7, 400)
    #ax.set_xticks([10, 30, 100, 300])
    #ax.set_xticklabels([10, 30, 100, 300])
    #plt.minorticks_off()
    
    # Display or save
    plt.tight_layout()
    #plt.show()
    plt.savefig("camel_sed_maxwellian_late.png", dpi=300)
    plt.close()


if __name__=="__main__":
    ultralong()
