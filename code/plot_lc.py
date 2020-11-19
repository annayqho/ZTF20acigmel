""" Plot radio/mm light curve """

from get_radio import *

# Font sizes
large = 14
medium = 11
small = 9


# Initializing plot
fig,ax = plt.subplots(1,1,figsize=(8,6))

# Get radio data
islim, tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()

# Plot the NOEMA LC
choose = freq==94.245
ax.errorbar(
        days[choose], flux[choose], 
        yerr=np.sqrt(eflux_form[choose]**2+eflux_sys[choose]**2), 
        fmt='o', c='k', label='NOEMA Band 1')

# Plot the ATCA 34 GHz LC
choose = np.logical_and(freq==34, islim==False)
ax.errorbar(
        days[choose], flux[choose], 
        yerr=np.sqrt(eflux_form[choose]**2+eflux_sys[choose]**2), 
        fmt='D', c='orange', label='ATCA 34 GHz')
choose = np.logical_and(freq==34, islim==True)
ax.scatter(days[choose][0], flux[choose][0], marker='_', c='orange')
ax.arrow(
        days[choose][0], flux[choose][0], 0, -flux[choose][0]/5, 
        length_includes_head=True, head_length=flux[choose][0]/10, 
        head_width=days[choose][0]/20, color='orange')

# Formatting
ax.set_xlabel("Days Since 2020 Oct 10", fontsize=large)
ax.set_ylabel("Flux Density", fontsize=large)
ax.set_ylim(0.07,1.2)
ax.set_yscale('log')
ax.tick_params(axis='both', labelsize=large)
ax.legend(loc='lower right', fontsize=medium)

# Display
plt.tight_layout()
plt.show()
