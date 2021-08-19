""" Fit a relativistic Maxwellian """

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset,inset_axes
from get_radio import *
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15
from format import *

d = get_format()


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


def at2020xnd(ax,col='#ef5675'):
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
    marker = '*'
    msize = 14
    ax.errorbar(x, y, ey, 
            fmt='%s' %marker, c=col, ms=msize, label=None)
    ax.errorbar(
            -100,-100,0,fmt='%s' %marker,c=col,
            ms=msize,label='AT2020xnd (40d)')

    # Fit for a Maxwellian w/o physical parameters
    p0 = [0.0003, 5.5E4, 0.7]

    # Fit for the Maxwellian in terms of physical quantities
    popt, pcov = curve_fit(fitfunc, x[1:], y[1:], p0=p0, #maxfev=1000000,
            sigma=ey[1:], absolute_sigma=True, 
            bounds=((0.0001,1E4,0.1),(0.0005,7E4,1.2)))
    xfit = np.logspace(0,2.5,200)
    yfit = fitfunc(xfit, *popt)
    print("AT2020xnd Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))

    #yfit = fitfunc(xfit, *p0)
    yfit = fitfunc(xfit, *popt)
    ax.plot(
            xfit,yfit, ls='-', zorder=5, color=col,
            label=r'Maxwellian ($\nu_m\approx1$ GHz)')

    # Plot a power law
    xplot = np.linspace(100,700)
    yplot = y[4]*(xplot/x[4])**(-1.5)
    ax.plot(xplot, yplot, c=col, ls='--', lw=1, label=None)

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
        choose = np.logical_and.reduce((days>b-b/20, days<b+b/20, islim==False))
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

        # Exclude points below 9 GHz
        keep = x >= 9
        x = x[keep]
        y = y[keep]
        ey = ey[keep]

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
        print("Maxwellian fit, Day %s:" %b)
        for i,param in enumerate(popt):
            print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))

        yfit = fitfunc(xfit, *popt)
        ax.plot(
                xfit,yfit, c='k', ls='--', zorder=5)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks([10,20,30,50,100])
        ax.set_xticklabels([10,20,30,50,100])
        ax.set_yticks([0.05,0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([0.05,0.2, 0.4, 0.6, 0.8])
        plt.minorticks_off()
        ax.set_xlim(10, 110)
        ax.set_ylim(0.05, 0.45)
        ax.legend(fontsize=d['font_small'])


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
    ax.plot(xfit,yfit, c='k', ls='-', zorder=5, label=None)



def at2018cow_panel(ax,factor,col='#7a5195'):
    """ factor: amount to scale the LC by """
    x = np.array([9.0, 34.0, 243.3, 259.3, 341.5, 357.5])
    y = np.array([0.27, 5.6, 36.6, 31.21, 19.49, 17.42])/factor
    ey = np.array([0.06, 0.16, 0.81, 0.92, 1.47, 2.8])/factor
    td = 10.5
    # Plot the data
    marker = 'o'
    msize = 8
    ax.errorbar(x, y, ey, 
            fmt='%s' %marker, c=col, label=None, ms=msize)
    ax.errorbar(-100, -100, 0, 
            fmt='%s' %marker, c=col, label='0.02xAT2018cow (10d)', ms=msize)

    # Fit for a Maxwellian w/o physical parameters
    p0 = [5E-2, 3E4, 1]
    popt, pcov = curve_fit(fitfunc, x, y, sigma=ey, absolute_sigma=True, p0=p0, maxfev=1000000)
    xfit = np.logspace(0,3,3000)
    yfit = fitfunc(xfit, *popt)
    print("Cow Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
    ax.plot(
            xfit,yfit, ls='-', zorder=5, 
            label=r'Maxwellian ($\nu_m\approx2$ GHz)', color=col)

    # Plot a power law
    xplot = np.linspace(250,1000)
    yplot = y[2]*(xplot/x[2])**(-1.5)
    ax.plot(xplot, yplot, c=col, ls='--', lw=1, label=None)


def at2018cow_plot():
    fig,ax = plt.subplots(1,1,figsize=(3.5,3))

    at2018cow_panel(ax,50)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xticks([10,20,30,50,100,200,300])
    ax.set_xticklabels([10,20,30,50,100,200,300])
    ax.set_yticks([0.2,0.4, 1, 3, 10, 30, 50])
    ax.set_yticklabels([0.2,0.4, 1, 3, 10, 30, 50])
    plt.minorticks_off()

    ax.set_xlabel(r"$\nu_\mathrm{rest}$ (GHz)", fontsize=d['font_med'])
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_med'])
    ax.tick_params(axis='both', labelsize=d['font_small'])

    ax.set_xlim(7, 400)
    ax.set_ylim(0.18, 80)

    plt.tight_layout()
    #plt.show()
    plt.savefig("cow_sed_maxwellian.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()


def ultralong_panel(ax, col='#ffa600'):
    z = 0.347
    x = np.array([4.8, 7.4, 9.5, 13.5, 16.0, 22.0]) * (1+z)
    y = np.array([237, 298, 293, 146, 89, 104]) / (1+z) / 1000
    ey = np.array([17, 12, 14, 20, 22, 13]) / (1+z) / 1000
    td = 2

    # Plot the data
    marker = 'D'
    msize = 8
    ax.errorbar(x, y, ey, 
            fmt='%s' %marker, c=col, label=None, ms=msize)
    ax.errorbar(-100, -100, 0, 
            fmt='%s' %marker, c=col, label='GRB130925A (2d)', ms=msize)

    # Fit for a Maxwellian w/o physical parameters
    p0 = [5E-1, 1000, 1E-1]

    popt, pcov = curve_fit(
            fitfunc, x[:-1], y[:-1], sigma=ey[:-1], absolute_sigma=True, 
            p0=p0, maxfev=100000)
    xfit = np.linspace(1,30)
    yfit = fitfunc(xfit, *popt)
    print("GRB Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
    ax.plot(
            xfit,yfit, ls='-', zorder=5, 
            label=r'Maxwellian ($\nu_m\approx0.1$ GHz)', c=col)

    xplot = np.linspace(12.5,40)
    yplot = y[2]*(xplot/x[2])**(-1.5)
    ax.plot(xplot, yplot, c=col, lw=1, ls='--', label=None)


def css161010_panel(ax, col='#ffa600'):
    factor = 6
    z = 0.034
    dat = np.loadtxt("css161010_radio_sed.txt",dtype=float,delimiter=',')
    x = dat[:,0]/1E9
    y = dat[:,1]*1E3
    print(x,y)

    # Plot the data
    marker = 'h'
    ax.scatter(x, y/factor, marker=marker , c=col, 
            label='0.2xCSS161010 (99d)', s=50)

    # Maxwellian fit
    p0 = [0.0003, 5.5E4, 0.7]
    popt, pcov = curve_fit(fitfunc, x, y, p0=p0)

    xfit = np.linspace(0.1,100,2000)
    yfit = fitfunc(xfit, *popt)
    print("CSS Maxwellian fit:")
    for i,param in enumerate(popt):
        print("%s +/- %s" %(param,np.sqrt(pcov[i,i])))
    ax.plot(
            xfit,yfit/factor, ls='-', zorder=5, 
            label=r'Maxwellian ($\nu_m\approx%s$ GHz)' %np.round(popt[2],1), 
            c=col)

    # Plot a power law as well
    xplot = np.linspace(10,100)
    yplot = y[20]*(xplot/x[20])**(-1.5)
    ax.plot(xplot, yplot/factor, c=col, lw=1, ls='--', label=None)


def ultralong():
    fig,ax = plt.subplots(1,1,figsize=(3.5,3))
    ultralong_panel(ax)
    ax.legend(fontsize=d['font_small'])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(10,250)
    ax.set_xlim(1,34)
    ax.set_yticks([50,70,100,150,200])
    ax.set_yticklabels([50,70,100,150,200])
    ax.set_xticks([5,7,10,15,20,30])
    ax.set_xticklabels([5,7,10,15,20,30])
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_med'])
    ax.tick_params(axis='both', labelsize=d['font_small'])
    ax.set_xlabel(r"$\nu_\mathrm{rest}$ (GHz)", fontsize=d['font_med'])
    plt.minorticks_off()
    plt.tight_layout()
    plt.show()
    #plt.savefig("ultralong_maxwellian.png", dpi=300, bbox_inches='tight',
    #        pad_inches=0.1)
    #plt.close()


def camel_late():
    """ Generate plots of the Camel """

    # one-panel
    fig,ax = plt.subplots(1,1,figsize=(3.5,3), sharey=True)

    at2020xnd_low_freq_late(ax)

    # Formatting
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_med'])
    ax.tick_params(axis='both', labelsize=d['font_small'])
    ax.set_xlabel(r"$\nu_{\mathrm{rest}}$ (GHz)", fontsize=d['font_med'])
    ax.set_yticks([0.1,0.2,0.3,0.4])
    ax.set_yticklabels([0.1,0.2,0.3,0.4])

    # Display or save
    plt.tight_layout()
    #plt.show()
    plt.savefig("camel_sed_maxwellian_late.png", dpi=300, bbox_inches='tight',
            pad_inches=0.1)
    plt.close()


def camel_early_full():
    """ Generate plots of the Camel """

    # one-panel
    fig,ax = plt.subplots(1,1,figsize=(3.5,2.5), sharey=True)

    at2020xnd(ax)

    # Formatting
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_med'])
    ax.tick_params(axis='both', labelsize=d['font_small'])
    #ax.set_ylim(7E-2, 1)
    #ax = axarr[1]
    ax.set_xlabel(r"$\nu_{\mathrm{rest}}$ (GHz)", fontsize=d['font_med'])
    ax.set_yticks([0.2,0.4,0.6,0.8])
    ax.set_yticklabels([0.2,0.4,0.6,0.8])
    
    # Display or save
    plt.tight_layout()
    #plt.show()
    plt.savefig("camel_sed_maxwellian_full.png", dpi=300, bbox_inches='tight',
            pad_inches=0.1)
    plt.close()


def camel_early():
    """ Generate plots of the Camel """

    # one-panel
    fig,ax = plt.subplots(1,1,figsize=(3.5,3), sharey=True)

    at2020xnd_high_freq(ax)

    # Formatting
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_med'])
    ax.tick_params(axis='both', labelsize=d['font_small'])
    ax.set_ylim(7E-2, 1.2)
    ax.set_xlabel(r"$\nu_{\mathrm{rest}}$ (GHz)", fontsize=d['font_med'])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_yticks([0.1,0.2,0.4,0.6,1.0])
    ax.set_yticklabels([0.1,0.2,0.4,0.6,1.0])
    ax.set_xticks([100,150,200,300])
    ax.set_xticklabels([100,150,200,300])
    ax.set_xlim(80,300)
    ax.legend(loc='upper right', fontsize=d['font_small'])
    
    # Display or save
    plt.tight_layout()
    plt.savefig("camel_sed_maxwellian_high_freq.png", dpi=200, bbox_inches='tight',
            pad_inches=0.1)
    #plt.show()
    plt.close()


def combined():
    """ combine panels into one big figure """
    fig,ax = plt.subplots(1,1,figsize=(6,5))

    # First, plot the AT2020xnd part
    at2020xnd(ax)

    # Then, plot the AT2018cow part
    at2018cow_panel(ax,50)

    # Then, plot the ultra-long GRB part
    ultralong_panel(ax,col='k')

    # Now CSS161010
    css161010_panel(ax)

    # Power law for legend
    ax.plot([0,1],[0,1], c='grey', ls='--', label=r'Power law ($\nu^{-1.5}$)', lw=1)

    # Formatting
    plt.legend(
        bbox_to_anchor=(0,1.02,1,1.02),loc='lower left',mode='expand',
        borderaxespad=0., ncol=2, fontsize=d['font_small'])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(4E-2, 1.5)
    ax.set_xticks([1,2,4,6,10,20,40,100,200,400])
    ax.set_xticklabels([1,2,4,6,10,20,40,100,200,400])
    ax.set_yticks([0.05,0.1,0.2,0.4,1])
    ax.set_yticklabels([0.05,0.1,0.2,0.4,1])
    ax.set_xlim(0.9, 600)
    ax.tick_params(axis='both', labelsize=d['font_med'])
    ax.set_xlabel(r"$\nu_{\mathrm{rest}}$ (GHz)", fontsize=d['font_large'])
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_large'])
    plt.tight_layout()

    # Save
    #plt.show()
    plt.savefig("maxwellian_combined.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    plt.close()

if __name__=="__main__":
    camel_late()
