import numpy as np

def get_survey_params(survey='SPT'):
    sensitivity = 0
    area = 0
    if survey=='SPT':
        # For SPT
        sensitivity = 15 # mJy
        area = 1500 / 41253 # as a fraction of the sky
    elif survey=='ACT':
        # For ACT
        sensitivity = 15*6 # mJy
        area = 0.4 # as a fraction of the sky
    elif survey=='S4wide':
        # For S4-wide (Collaboration et al. 2019), one week
        sensitivity = 3*6 # mJy
        area = 0.5 # as a fraction of the sky
    elif survey=='S4deep':
        # For S4-deep
        sensitivity = 0.8*6 # mJy, from Collaboration et al. 2019
        area = 0.03 # as a fraction of the sky
    elif survey=='Simons':
        # For Simons
        sensitivity = 4*6 # mJy, from Collaboration et al. 2019
        area = 0.4 # as a fraction of the sky
    return sensitivity, area

def get_source_params(source='AT2020xnd'):
    lpeak = 0
    rate = 0
    if source=='AT2020xnd':
        # AT2020xnd parameters
        lpeak = 2E30
        rate = 0.001*7E-5 # /Mpc3/yr
    elif source=='LGRB':
        # LGRB parameters
        lpeak = 1E32
        rate = 0.25 * 4.2E-10
    elif source=='LLGRB':
        # LLGRB parameters
        lpeak = 1E29
        rate = 2.3E-7
    elif source=='SN':
        # SN parameters
        lpeak = 1E27
        rate = 7E-5
    return lpeak,rate

print("Calculation for SN")
surveys = ['SPT', 'ACT', 'S4wide', 'S4deep']
lpeak,rate = get_source_params(source='AT2020xnd')
for s in surveys:
    print(s)

    ### Get parameters
    sensitivity, area = get_survey_params(survey=s)

    ### Calculate volume
    sensitivity_cgs = sensitivity * 1E-3 * 1E-23
    d_lim = np.sqrt(lpeak / (4*np.pi*sensitivity_cgs))
    d_lim_mpc = d_lim / 3.086E24
    print("Volume: ", d_lim_mpc)

    ### Calculate rate in that volume
    det_rate = (4/3) * np.pi * d_lim_mpc**3 * rate * area
    print("Rate: ", det_rate)
