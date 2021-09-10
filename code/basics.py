import numpy as np
from astropy.cosmology import Planck15
from get_radio import *


def get_z():
    return 0.2433


def get_dL():
    """ This number is reported at the end of Section 1 in the paper """
    z = get_z()
    dL = Planck15.luminosity_distance(z=z).value
    return dL
    

def get_dcm():
    """ This number is reported at the end of Section 1 in the paper """
    z = get_z()
    dcm = Planck15.luminosity_distance(z=z).cgs.value
    return dcm


def get_mm_peak():
    """ This number is reported in Section 2 """
    islim, tel, freq, days, flux, eflux = get_data_all()
    choose = np.logical_and(freq>70, freq<80)
    ind = np.argmax(flux[choose])
    peak_flux = flux[choose][ind]
    peak_eflux = eflux[choose][ind]
    print(peak_flux, peak_eflux)
    
    dcm = get_dcm()
    peak_lum = peak_flux * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    peak_elum = peak_eflux * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    print(peak_lum/1E30, peak_elum/1E30)


def get_cm_peak():
    """ This number is reported in Section 2 """
    islim, tel, freq, days, flux, eflux = get_data_all()
    choose = np.logical_and.reduce((freq>9, freq<11, islim==False))
    ind = np.argmax(flux[choose])
    peak_flux = flux[choose][ind]
    peak_eflux = eflux[choose][ind]
    print(peak_flux, peak_eflux)
    
    dcm = get_dcm()
    peak_lum = peak_flux * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    peak_elum = peak_eflux * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    print(peak_lum/1E29, peak_elum/1E29)


if __name__=="__main__":
    #get_dL()
    #get_mm_peak()
    get_cm_peak()
