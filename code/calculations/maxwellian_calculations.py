import numpy as np

# Note: the equations are only correct for D=1266 Mpc.
# So, you need the flux densities to be corrected to that luminosity distance.

# AT2020xnd values
fm = 0.00025978580596915037 # mJy
taum = 61690.7749196531
nuT = 0.7202302548516042 # GHz
t = 40

# AT2018cow values
d_mpc = 60
fm = 0.03813173754390941  * (d_mpc/1261)**2 # equations are for 20xnd
taum = 19486.806270613346
nuT = 1.8801301063938891
t = 10

# CSS161010 values
d_mpc = 154
fm = 0.031638890979651516 * (d_mpc/1261)**2
taum = 725.4236839560893
nuT = 0.15988247897298238
t = 99
# 
# GRB 130925A values
d_mpc = 1895
fm = 0.15825289310156195 * (d_mpc/1266)**2
taum = 17360.69575983693
nuT = 0.13005043044403894
t = 2
# 
# # AT2020xnd Maxwellian fit, Day 71:
# fm = 0.0005999999999999996 
# taum = 4476.715522437825 
# nuT = 0.5574496313370073 
# t = 71/1.2442
# 
# # Maxwellian fit, Day 95:
# fm = 0.0005999999999999998 
# taum = 2189.284336888685 
# nuT = 0.5582891494084399 
# t = 95/1.2442
# 
# # Maxwellian fit, Day 132:
# fm = 0.00016644499735265936 
# taum = 3971.713887211959 
# nuT = 0.33598250530435236
# t = 132/1.2442

# Begin Calculation

nuT_MHz = nuT*1000
t50 = t/50

v = ((fm/nuT_MHz**2)*(8.3**2/7.45E-10)*(1/t50**2))**(1/4)
print(v)
print("Velocity: %sc" %np.round(v*0.1,1))

B = (nuT_MHz/8.30)/v**4
print("Magnetic field strength: %s G" %np.round(B,10))

ne = 100*(taum/3.73E7)*B*v**9/t50
print("Ambient density: %s cm-3" %np.round(ne))
