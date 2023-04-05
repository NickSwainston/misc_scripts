import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

from pulsar_spectra.models import low_frequency_turn_over_power_law

# freq range
freqs_MHz = np.logspace(np.log10(100), np.log10(2200), num=100)

# make dummy lfto function
flux = low_frequency_turn_over_power_law(freqs_MHz, 180, -2, 1, 2.1, 500)

# Set up default mpl markers
plotsize = 3.2

# Set up plot
fig, ax = plt.subplots(figsize=(plotsize, plotsize))
# plot data
ax.plot(freqs_MHz, flux, 'black', linewidth=3, zorder=0.5)

# Format plot and save
ax.set_xscale('log')
ax.set_yscale('log')
ax.get_xaxis().set_major_formatter(FormatStrFormatter('%g'))
ax.get_yaxis().set_major_formatter(FormatStrFormatter('%g'))
ax.tick_params(which='both', direction='in', top=1, right=1)
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Flux Density (Arbitrary Units)')

plt.savefig("lfto_function.png", bbox_inches='tight', dpi=300)
plt.clf()
