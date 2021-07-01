""" Plot spectral index between various bands as a function of time """
import matplotlib.pyplot as plt
import numpy as np

fig,ax = plt.subplots(1,1,figsize=(6,4.5))

# Highest: NOEMA to SMA (146 to 185)
# c = '#ffa600'
# dt = 46.5
# nu1 = 146
# nu2 = 185
# f1 = 0.153
# f2 = 0.150
# alpha = np.log(f1/f2)/np.log(nu1/nu2)
# plt.scatter(dt, alpha, c=c, marker='s', label='146/185')
# plt.arrow(dt, alpha, 0, -1, color=c, length_includes_head=True,
#         head_length=0.5, head_width=3)
# 
# Highest: NOEMA to SMA (94 to 185)
#dt = 22.5
#nu1 = 94
#nu2 = 185
#f1 = 0.648
#f2 = 0.48
#alpha = np.log(f1/f2)/np.log(nu1/nu2)
#plt.scatter(dt, alpha, edgecolor=c, facecolor='white', marker='s', label='94/185')
#plt.arrow(dt, alpha, 0, -1, color=c, length_includes_head=True,
#        head_length=0.5, head_width=3)

# NOEMA 131 to 146
c = 'k'
nu1 = 131
nu2 = 146
dt = 46
f1 = 0.318
f2 = 0.179
ef1 = 0.049
ef2 = 0.056
alpha = np.log(f1/f2)/np.log(nu1/nu2)
ealpha = (1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)
#plt.errorbar(dt, alpha, yerr=ealpha, fmt='o', c=c, mfc='white', mec=c, label='131/146', lw=0.5)

# NOEMA 94 to 131
c = 'k'
nu1 = 94
nu2 = 131
dt = 46
f1 = 0.569
f2 = 0.318
ef1 = 0.042
ef2 = 0.049
alpha = np.log(f1/f2)/np.log(nu1/nu2)
ealpha = (1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)
plt.errorbar(
        dt, alpha, yerr=ealpha, fmt='o', c=c, 
        mfc='white', mec=c, label='94/131', lw=0.5)

# NOEMA 79 to 94
c = '#ffa600'
nu1 = 79
nu2 = 94

dt = [17]
f1 = 0.392
f2 = 0.305
ef1 = 0.059
ef2 = 0.057
alpha = [np.log(f1/f2)/np.log(nu1/nu2)]
ealpha = [(1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)]

dt.append(24)
f1 = 0.679
f2 = 0.648
ef1 = 0.046
ef2 = 0.044
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

dt.append(31)
f1 = 1.09
f2 = 1.03
ef1 = 0.049
ef2 = 0.044
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

dt.append(39)
f1 = 0.924
f2 = 0.868
ef1 = 0.050
ef2 = 0.046
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

dt.append(46)
f1 = 0.822
f2 = 0.569
ef1 = 0.048
ef2 = 0.042
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

dt.append(67)
f1 = 0.247
f2 = 0.130
ef1 = 0.037
ef2 = 0.036
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

ax.scatter(dt, alpha, marker='o', c=c, label="79/94")
ax.errorbar(dt, alpha, yerr=ealpha, fmt='o-', c=c, lw=0.5)

# NOEMA 79 to VLA 36.2
c = '#ff6e54'
nu1 = 36.2
nu2 = 79
dt = [37.5]
f1 = 0.675
f2 = 0.924
ef1 = 0.171
ef2 = 0.050
alpha = [np.log(f1/f2)/np.log(nu1/nu2)]
ealpha = [(1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)]
dt.append(48.5)
f1 = 0.668
f2 = 0.822
ef1 = 0.168
ef2 = 0.048
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

ax.scatter(dt, alpha, marker='s', c=c, label="36.2/79")
ax.errorbar(dt, alpha, yerr=ealpha, fmt='s-', c=c, lw=0.5)

# VLA/ATCA 27 to VLA 36.2
c = '#dd5182'
nu1 = 27
nu2 = 36.2
dt = [36]
f1 = 0.497
f2 = 0.675
ef1 = 0.011
ef2 = 0.171
alpha = [np.log(f1/f2)/np.log(nu1/nu2)]
ealpha = [(1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)]
dt.append(51.9)
f1 = 0.621
f2 = 0.668
ef1 = 0.132
ef2 = 0.168
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))
dt.append(71.9)
f1 = 0.450
f2 = 0.209
ef1 = 0.096
ef2 = 0.063
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

ax.scatter(
        dt, alpha, marker='P', c=c, label="33/45", s=50)
ax.errorbar(dt, alpha, yerr=ealpha, fmt='P-', c=c, lw=0.5, ms=8)

# VLA 17.7 to VLA/ATCA 27
c = '#955196'
nu1 = 17.7
nu2 = 27
dt = [71]
f1 = 0.484
f2 = 0.450
ef1 = 0.010
ef2 = 0.096
alpha = [np.log(f1/f2)/np.log(nu1/nu2)]
ealpha = [(1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)]
dt.append(94)
f1 = 0.301
f2 = 0.213
ef1 = 0.065
ef2 = 0.048
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))
dt.append(131)
f1 = 0.087
f2 = 0.048
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append(0)
ax.scatter(dt, alpha, marker='X', c=c, label="22/33")
ax.errorbar(dt, alpha, yerr=ealpha, fmt='X-', c=c, lw=0.5, ms=8)

# 12 GHz to 17.7 GHz
c = '#444e86'
nu1 = 12
nu2 = 17.7
dt = [71]
f1 = 0.401
f2 = 0.484
ef1 = 0.046
ef2 = 0.010
alpha = [np.log(f1/f2)/np.log(nu1/nu2)]
ealpha = [(1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1)]
dt.append(94)
f1 = 0.278
f2 = 0.301
ef1 = 0.032
ef2 = 0.065
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))
dt.append(131)
f1 = 0.122
f2 = 0.087
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))
ax.scatter(dt, alpha, marker='D', c=c, label="15/22")
ax.errorbar(dt, alpha, yerr=ealpha, fmt='D-', c=c, lw=0.5)


# VLA 8 to 12 GHz
c = '#003f5c'
nu1 = 8
nu2 = 12

dt = [17]
f1 = 0.027
f2 = 0.037
alpha = [np.log(f1/f2)/np.log(nu1/nu2)]
ealpha = [0]
ax.arrow(
        dt[0], alpha[0], 0, alpha[0]/3, length_includes_head=True, color=c,
        head_width=dt[0]/15, head_length=alpha[0]/5)

dt.append(25)
f1 = 0.057
f2 = 0.095
ef1 = 0.005
ef2 = 0.011
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

dt.append(71)
f1 = 0.180
f2 = 0.401
ef1 = 0.023
ef2 = 0.046
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

dt.append(94)
f1 = 0.168
f2 = 0.278
ef1 = 0.022
ef2 = 0.032
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))

dt.append(131)
f1 = 0.109
f2 = 0.122
ef1 = 0.010
ef2 = 0.016
alpha.append(np.log(f1/f2)/np.log(nu1/nu2))
ealpha.append((1/np.log(nu1/nu2)) * (1/(f1*f2)) * (f1*ef2-f2*ef1))
ax.scatter(dt, alpha, c=c, marker='*', label="10/15", s=50)
ax.errorbar(dt, alpha, yerr=ealpha, marker='*', c=c, lw=0.5, ms=10)

ax.set_xscale('log')
ax.set_ylim(-4.0,2.1)
ax.set_xlabel("Days since 2020 Oct 10.0", fontsize=16)
ax.set_ylabel(r"Spectral index $\alpha$ ($f_\nu \propto \nu^{\alpha}$)", fontsize=16)
plt.legend(loc='lower left', ncol=2, fontsize=10.5)
ax.tick_params(axis='both', labelsize=14)
ax.set_xticks([20,30,40,60,100])
ax.set_xticklabels([20,30,40,60,100])
plt.tight_layout()
#plt.show()
plt.savefig("spindex_time.png", dpi=200)
plt.close()

