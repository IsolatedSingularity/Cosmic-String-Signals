#%% Importing Modules
import numpy as np
import matplotlib.pyplot as plt
import powerbox as pbox

# %% 2D Gaussian field w/ power-law power-spectrum
pb = pbox.PowerBox(
    N=512,                     # Number of grid-points in the box
    dim=2,                     # 2D box
    pk = lambda k: 0.1*k**-2., # The power-spectrum
    boxlength = 1.0,           # Size of the box (sets the units of k in pk)
    seed = 1010,               # Set a seed to ensure the box looks the same every time (optional)
    ensure_physical=True       # ** Ensure the delta_x is a physically valid over-density **
)


plt.imshow(pb.delta_x(),extent=(0,1,0,1))
plt.colorbar()
plt.show()

# %% 2D log-normal field power-law power-spectrum
lnpb = pbox.LogNormalPowerBox(
    N=512,                     # Number of grid-points in the box
    dim=2,                     # 2D box
    pk = lambda k: 0.1*k**-2., # The power-spectrum
    boxlength = 1.0,           # Size of the box (sets the units of k in pk)
    seed = 1010                # Use the same seed as our powerbox
)
plt.imshow(lnpb.delta_x(),extent=(0,1,0,1))
plt.colorbar()
plt.show()

# %% Discrete field sample

fig, ax = plt.subplots(1,2, sharex=True,sharey=True,gridspec_kw={"hspace":0}, subplot_kw={"ylim":(0,1),"xlim":(0,1)}, figsize=(10,5))

# Create a discrete sample using the PowerBox instance.
samples = pb.create_discrete_sample(nbar=50000,      # nbar specifies the number density
                                    min_at_zero=True  # by default the samples are centred at 0. This shifts them to be positive.
                                   )
ln_samples = lnpb.create_discrete_sample(nbar=50000, min_at_zero=True)

# Plot the samples
ax[0].scatter(samples[:,0],samples[:,1], alpha=0.2,s=1)
ax[1].scatter(ln_samples[:,0],ln_samples[:,1],alpha=0.2,s=1)
plt.show()


# %% Plotting power spectrum
from powerbox import get_power

# Only two arguments required when passing a field
p_k_field, bins_field = get_power(pb.delta_x(), pb.boxlength)
p_k_lnfield, bins_lnfield = get_power(lnpb.delta_x(), lnpb.boxlength)

# The number of grid points are also required when passing the samples
p_k_samples, bins_samples = get_power(samples, pb.boxlength,N=pb.N)
p_k_lnsamples, bins_lnsamples = get_power(ln_samples, lnpb.boxlength,N=lnpb.N)


plt.plot(bins_field, 0.1*bins_field**-2., label="Input Power")

plt.plot(bins_field, p_k_field,label="Normal Field Power")
plt.plot(bins_samples, p_k_samples,label="Normal Sample Power")
plt.plot(bins_lnfield, p_k_lnfield,label="Log-Normal Field Power")
plt.plot(bins_lnsamples, p_k_lnsamples,label="Log-Normal Sample Power")

plt.legend()
plt.xscale('log')
plt.yscale('log')



# %%
