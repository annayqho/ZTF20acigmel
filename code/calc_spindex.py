""" Calculate various spectral indices """

import numpy as np

# The optically thin spectral index on Day 21
f1 = 0.234
f2 = 0.146
ef1 = 0.034
ef2 = 0.033
nu1 = 79
nu2 = 94

alpha = np.log(f1/f2)/np.log(nu1/nu2)
ealpha = (1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)

print(alpha, ealpha)

# The optically thin spectral index on Day 28
f1 = 0.91
f2 = 0.87
ef1 = 0.04
ef2 = 0.03
nu1 = 79
nu2 = 94

alpha = np.log(f1/f2)/np.log(nu1/nu2)
ealpha = (1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)

print(alpha, ealpha)

# The optically thin index when you have the SMA upper limit
f1 = 0.57
f2 = 0.95
nu1 = 185
nu2 = 94
alpha = np.log(f1/f2)/np.log(nu1/nu2)
print(alpha)
