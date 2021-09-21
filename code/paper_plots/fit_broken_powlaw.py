""" Fit a broken power law to the SED 

and determine the evolution of that SED with time
"""

from scipy.optimize import curve_fit
import sys
sys.path.append("/Users/annaho/Dropbox/astronomy/papers_active/ZTF20acigmel/code")
from get_radio import *
from format import *
from basics import get_z

d = get_format()

z = get_z()
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
    beta1 = 2.5
    beta2 = -1
    nu,t = x_in

    pref = Fa0 * (t/t0)**(alpha1)
    bot = nua0 * (t/t0)**(alpha2)
    Fnu = pref * ((nu/bot)**(-s*beta1) + (nu/bot)**(-s*beta2))**(-1/s)

    return Fnu


def func_no_spindex_const_shock(x_in, Fa0, nua0, alpha1, s):
    # for a constant velocity, target_value = 1
    # for a decelerating shock, target_value = 0.8
    target_value = 1
    alpha2 = 9*alpha1/19-target_value
    beta1 = 2.5
    beta2 = -1
    nu,t = x_in

    pref = Fa0 * (t/t0)**(alpha1)
    bot = nua0 * (t/t0)**(alpha2)
    Fnu = pref * ((nu/bot)**(-s*beta1) + (nu/bot)**(-s*beta2))**(-1/s)

    return Fnu


def func_spec(nu, Fa0, nua0, beta1, beta2, s):
    Fnu = Fa0 * ((nu/nua0)**(-s*beta1) + (nu/nua0)**(-s*beta2))**(-1/s)
    return Fnu


def func_spec_no_spindex(nu, Fa0, nua0):
    s = 1
    beta1 = 5/2
    beta2 = -1
    Fnu = Fa0 * ((nu/nua0)**(-s*beta1) + (nu/nua0)**(-s*beta2))**(-1/s)
    #Fnu = Fa0 * ((nu/nua0)**(beta1) + (nu/nua0)**(beta2))
    return Fnu


def run_late_time():
    """ Run the full process for the late-time SEDs, with
    power-law evolution """
    islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()

    # K-correct and convert to rest-frame before modeling
    freq = freq_obs[islim==False] * (1+z)
    flux = flux_obs[islim==False] / (1+z)
    eflux = eflux_obs[islim==False] / (1+z)
    days = days_obs[islim==False] / (1+z)

    # Times of observation
    bins_obs = np.array([13, 18, 24, 28.3, 31.8, 38, 46, 51.9, 71, 95, 132])
    bins = bins_obs / (1+z)

    # Only plot the last three epochs
    use_ind = np.arange(len(bins))[8:]

    t = [] # time 
    x = [] # frequency
    y = [] # flux
    ey = []
    for ii in use_ind:
        bin = bins[ii]
        choose = np.abs(days-bin) < bin/20
        [t.append(bin) for i in np.arange(sum(choose))]
        [x.append(val) for val in freq[choose]]
        [y.append(val) for val in flux[choose]]
        [ey.append(val) for val in eflux[choose]]
    t = np.array(t)
    x = np.array(x)
    y = np.array(y)
    ey = np.array(ey)

    # Plot the data
    fig,axarr = plt.subplots(1,2,figsize=(6,2.5),sharex=True, sharey=True)
    col = d['colors']['9']
    col = col[::-1][6:]
    markers = ['o', 's', 'D']
    msize = [6, 6, 6]

    for ax in axarr:
        j = 0
        for ii in use_ind:
            lab = int(bins_obs[ii])
            choose = np.logical_and(t==bins[ii], x>6*(1+z))
            ax.errorbar(
                    x[choose], y[choose], yerr=ey[choose], 
                    fmt=markers[j], c=col[j], ms=msize[j],
                    label='%s d' %lab)
            choose = np.logical_and(t==bins[ii], x==6*(1+z))
            ax.errorbar(
                    x[choose], y[choose], yerr=ey[choose], 
                    fmt='o', mfc='white', mec=col[j], c=col[j])
            j += 1

    # For the fitting, ignore the 6 GHz point
    tofit = x > (6*(1+z)) # rest-frame
    x = x[tofit]
    y = y[tofit]
    t = t[tofit]
    ey = ey[tofit]

    # The first fit does not assume a constant shock.
    print("Doing the first fit: no assumptions")
    ax = axarr[0]

    p0 = np.array([0.86, 17, -2.2, -0.9, 1])
    popt, pcov = curve_fit(
            func_no_spindex, (x, t), y, 
            p0=p0, sigma=ey, absolute_sigma=True, maxfev=10000)
            #bounds=([0.1,10,-4,-4,0.1,-4],[2,50,-0.1,-0.1,4,-0.1]))
    for i in np.arange(len(p0)):
        print("%s +/- %s" %(popt[i], np.sqrt(pcov[i,i])))

    nuplot = np.logspace(0.5,2.4)
    j = 0
    for ii in use_ind:
        lab = int(bins_obs[ii]/(1+z))
        tplot = np.array([bins[ii]]*len(nuplot))
        fplot = func_no_spindex((nuplot,tplot), *popt)
        print(bins_obs[ii]/(1+z))
        ax.plot(nuplot, fplot, c=col[j])
        j += 1

    # Calculate the chi squared
    chisq = 0
    for ii in use_ind:
        choose = t==bins[ii]
        f_model = func_no_spindex((x[choose],t[choose]), *popt)
        chisq += sum(((f_model - y[choose])/ey[choose])**2)
    print("number of DOF:")
    print(len(t)-len(p0))
    red_chisq = np.round(chisq / (len(t)-len(p0)),1)
    ax.text(0.05,0.9,r'$\chi^2_r\approx%s$' %red_chisq, 
            fontsize=d['font_small'], transform=ax.transAxes)


    # The second fit does assume a constant shock
    print("Doing the second fit: assume constant shock")
    ax = axarr[1]
    #p0 = np.array([0.55, 30, -7, 1]) #if you fit for just k
    p0 = np.array([0.86, 17, -0.9, 1])

    popt, pcov = curve_fit(
            func_no_spindex_const_shock, (x, t), y, 
            p0=p0, sigma=ey, absolute_sigma=True, maxfev=10000)
            #bounds=([0.1,10,-4,-4,0.1,-4],[2,50,-0.1,-0.1,4,-0.1]))
    for i in np.arange(len(p0)):
        print("%s +/- %s" %(popt[i], np.sqrt(pcov[i,i])))

    nuplot = np.logspace(0.5,2.4)
    j = 0
    for ii in use_ind:
        lab = int(bins_obs[ii]/(1+z))
        tplot = np.array([bins[ii]]*len(nuplot))
        fplot = func_no_spindex_const_shock((nuplot,tplot), *popt)
        ax.plot(nuplot, fplot, c=col[j])
        j += 1

    # Calculate the chi squared
    chisq = 0
    for ii in use_ind:
        choose = t==bins[ii]
        f_model = func_no_spindex_const_shock((x[choose],t[choose]), *popt)
        chisq += sum(((f_model - y[choose])/ey[choose])**2)
    print("number of DOF:")
    print(len(t)-len(p0))
    red_chisq = np.round(chisq / (len(t)-len(p0)),1)
    ax.text(0.05,0.9,r'Constant shock' %red_chisq, 
            fontsize=d['font_small'], transform=ax.transAxes)
    ax.text(0.05,0.8,r'$\chi^2_r\approx%s$' %red_chisq, 
            fontsize=d['font_small'], transform=ax.transAxes)

    for ax in axarr:
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(6.2, 124)
        ax.set_ylim(5E-2, 0.6)
        ax.set_xlabel(r"$\nu_{\mathrm{rest}}$ (GHz)", fontsize=d['font_med'])
    handles, labels = axarr[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3)
    plt.subplots_adjust(wspace=0.1)
    axarr[0].set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_med'])
    #plt.tight_layout()
    #plt.show()
    plt.savefig("broken_powlaw_fits.png", dpi=200, bbox_inches='tight',
            pad_inches=0.1)
    plt.close()


def run_one_epoch():
    """ run for just one epoch """
    islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()

    # Put everything into the rest-frame before modeling
    freq = freq_obs[islim==False] * (1+z)
    flux = flux_obs[islim==False] / (1+z)
    eflux = eflux_obs[islim==False] / (1+z)
    days = days_obs[islim==False] / (1+z)

    # Times of observation
    bins_obs = np.array([13, 18, 24, 28.3, 31.8, 38, 46, 51.9, 71, 95, 132])
    bins = bins_obs / (1+z)

    # Only plot the last three epochs
    use_ind = np.arange(len(bins))[6:7]

    t = [] # time 
    x = [] # frequency
    y = [] # flux
    ey = []
    for ii in use_ind:
        bin = bins[ii]
        choose = np.abs(days-bin) < bin/20
        [t.append(bin) for i in np.arange(sum(choose))]
        [x.append(val) for val in freq[choose]]
        [y.append(val) for val in flux[choose]]
        [ey.append(val) for val in eflux[choose]]
    t = np.array(t)
    x = np.array(x)
    y = np.array(y)
    ey = np.array(ey)

    # Plot the data
    fig,ax = plt.subplots(1,1,figsize=(3.5,2.5),sharex=True, sharey=True)
    col = d['colors']['3']
    j = 0
    for ii in use_ind:
        choose = np.logical_and(t==bins[ii], x>6*(1+z))
        ax.errorbar(
                x[choose], y[choose], yerr=ey[choose], fmt='o', c=col[j])
        choose = np.logical_and(t==bins[ii], x==6*(1+z))
        ax.errorbar(
                x[choose], y[choose], yerr=ey[choose], 
                fmt='o', mfc='white', mec=col[j], c=col[j])
        j += 1

    # For the fitting, ignore the 6 GHz point
    tofit = x > (6*(1+z)) # rest-frame
    x = x[tofit]
    y = y[tofit]
    t = t[tofit]
    ey = ey[tofit]

    # The first fit does not assume a constant shock.
    p0 = np.array([0.86, 17, 1])
    p0 = np.array([1, 17])
    popt, pcov = curve_fit(
            func_spec_no_spindex, x, y, 
            p0=p0, sigma=ey, absolute_sigma=True, maxfev=10000)
            #bounds=([0.1,10,-4,-4,0.1,-4],[2,50,-0.1,-0.1,4,-0.1]))
    for i in np.arange(len(p0)):
        print("%s +/- %s" %(popt[i], np.sqrt(pcov[i,i])))

    nuplot = np.logspace(0.5,2.4)
    j = 0
    for ii in use_ind:
        lab = int(bins_obs[ii]/(1+z))
        tplot = np.array([bins[ii]]*len(nuplot))
        fplot = func_spec_no_spindex(nuplot, *popt)
        print(bins_obs[ii]/(1+z))
        ax.plot(nuplot, fplot, c=col[j], label="%s d" %lab)
        j += 1

    # Calculate the chi squared
    chisq = 0
    for ii in use_ind:
        choose = t==bins[ii]
        f_model = func_spec_no_spindex(x[choose], *popt)
        chisq += sum(((f_model - y[choose])/ey[choose])**2)
    red_chisq = np.round(chisq / (len(t)-len(p0)),1)
    ax.text(0.05,0.7,r'$\chi^2_r\approx%s$' %red_chisq, 
            fontsize=d['font_small'], transform=ax.transAxes)


    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend(loc='upper left', ncol=3, fontsize=d['font_small'])
    ax.set_xlim(6.2, 124)
    ax.set_ylim(0.02, 0.6)
    ax.set_xlabel(r"$\nu_{\mathrm{rest}}$ (GHz)", fontsize=d['font_med'])
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_med'])
    plt.tight_layout()
    plt.show()



def run_early_time():
    """ Run the full process for the early-time SEDs, with
    power-law evolution """
    islim, tel, freq_obs, days_obs, flux_obs, eflux_obs = get_data_all()

    # Put everything into the rest-frame before modeling
    freq = freq_obs[islim==False] * (1+z)
    flux = flux_obs[islim==False] / (1+z)
    eflux = eflux_obs[islim==False] / (1+z)
    days = days_obs[islim==False] / (1+z)

    # Times of observation
    bins_obs = np.array([18, 24, 30.8, 38, 46, 51.9, 71, 95, 132])
    bins = bins_obs / (1+z)

    # Only plot the last three epochs
    use_ind = np.arange(len(bins))[0:6]

    # Plot the data
    fig,ax = plt.subplots(1,1,figsize=(3.5,3),sharex=True, sharey=True)
    col = d['colors']['9'][::-1]

    chisq = 0
    npt = 0
    for j,ii in enumerate(use_ind):
        bin = bins[ii]
        print(bin)
        choose = np.abs(days-bin) < bin/20
        x = freq[choose]
        y = flux[choose]
        ey = eflux[choose]

        choose = x>6*(1+z)
        ax.errorbar(
                x[choose], y[choose], yerr=ey[choose], fmt='o', c=col[j])
        choose = x==6*(1+z)
        ax.errorbar(
                x[choose], y[choose], yerr=ey[choose], 
                fmt='o', mfc='white', mec=col[j], c=col[j])

        # For the fitting, ignore the 6 GHz point
        tofit = x > (6*(1+z)) # rest-frame
        x = x[tofit]
        y = y[tofit]
        ey = ey[tofit]
        npt += len(x)

        print("Doing the fit")
        p0 = np.array([0.86, 17])
        popt, pcov = curve_fit(
                func_spec_no_spindex, x, y, 
                p0=p0, sigma=ey, absolute_sigma=True, maxfev=10000)
        for i in np.arange(len(p0)):
            print("%s +/- %s" %(popt[i], np.sqrt(pcov[i,i])))

        nuplot = np.logspace(0.5,2.4)
        tplot = np.array([bins[ii]]*len(nuplot))
        fplot = func_spec_no_spindex(nuplot, *popt)
        ax.plot(nuplot, fplot, c=col[j], label="%sd" %(int(bins_obs[ii])))

        # Calculate the chi squared
        f_model = func_spec_no_spindex(x, *popt)
        chisq += sum(((f_model - y)/ey)**2)
    red_chisq = np.round(chisq / (npt-2),1)
    ax.text(0.05,0.7,r'$\chi^2_r\approx%s$' %red_chisq, 
            fontsize=d['font_small'], transform=ax.transAxes)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend(loc='lower right', ncol=2, fontsize=d['font_small'])
    ax.set_xlim(6.2, 300)
    ax.set_ylim(1E-2, 0.6)
    ax.set_xlabel(r"$\nu_{\mathrm{rest}}$ (GHz)", fontsize=d['font_small'])
    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=d['font_small'])
    ax.tick_params(axis='both', labelsize=d['font_small'])
    plt.tight_layout()
    plt.show()
    #plt.savefig(
    #        "powlaw_fits_early.png", dpi=200, 
    #        bbox_inches='tight', pad_inches=0.1)
    #plt.close()


if __name__=="__main__":
    run_late_time()
