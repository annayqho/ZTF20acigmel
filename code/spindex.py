import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

def run(nu, f, ef):
    npts = len(nu)
    ndraw = 1000

    # Make 1000 new versions of a new flux array, repopulated with the flux
    # from the Gaussian dist with STD of that uncertainty

    f_new = np.zeros((ndraw,npts))
    for ii,val in enumerate(f):
        # log uncertainty is dx/x
        f_new[:,ii] = np.random.normal(loc=val, scale=ef[ii], size=ndraw)
    
    for val in f_new:
        plt.plot(nu,val,c='grey',alpha=0.2)

    alphas = np.zeros(ndraw)

    for ii in np.arange(ndraw):
        alphas[ii] = np.log10(f_new[ii][1]/f_new[ii][0])/np.log10(nu[1]/nu[0])

    # Measure the mean and standard deviation
    alpha = np.mean(alphas)
    ealpha = np.std(alphas)
