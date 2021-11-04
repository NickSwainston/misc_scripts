
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import math

from vcstools.metadb_utils import calc_ta_fwhm

def scale_for_all_freqs(ref_size, ref_freq, freq, a=-2.):
    # Scale these size estimates to all frequencies
    return freq**a * ref_freq**(-a) * ref_size

freq = np.logspace(math.log10(80) ,math.log10(300) , num=30)
#print(freq)

# Plot Full Radius Half Measures for different array phases
ctf = np.vectorize(calc_ta_fwhm)
compact_frhm = ctf(freq, array_phase='P2C') / 2 * 60 * 60
extended_frhm = ctf(freq, array_phase='P2E') / 2 * 60 * 60

legend_lines = []
#fig, ax1 = plt.subplots()
size = 5
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(2*size,size))
#ax1.plot(freq, compact_frhm, label="P2C FRHM")
legend_lines.append(ax1.plot(freq, extended_frhm, label="P2E FRHM", linestyle='dashed', color='purple'))
# fake line
legend_lines.append(ax1.plot(80, 60,color='white'))

sfaf = np.vectorize(scale_for_all_freqs)

# Chris's estimate median
chris_median_offset = 0.1304 * 60 #arc second
chris_median_freq = 200 # MHz
chris_median_scaled = sfaf(chris_median_offset, chris_median_freq, freq)
legend_lines.append(ax1.plot(freq, chris_median_scaled, label="Residual angular offset (50%)", linestyle='dotted', color='green'))
chris_median_response = np.exp(-(chris_median_scaled / extended_frhm)**2/(2))
ax2.plot(freq, chris_median_response, label="Residual angular offset (50%)", linestyle='dotted', color='green')

# Chris's esitmate 90%
chris_90_offset = 0.29105 * 60 #arc second
chris_90_freq = 200 # MHz
chris_90_scaled = sfaf(chris_90_offset, chris_90_freq, freq)
legend_lines.append(ax1.plot(freq, chris_90_scaled, label="Residual angular offset (90%)", linestyle='dotted', color='blue'))
chris_90_response = np.exp(-(chris_90_scaled / extended_frhm)**2/(2))
ax2.plot(freq, chris_90_response, label="Residual angular offset (90%)", linestyle='dotted', color='blue')

#plt.errorbar(freq, chris_scaled, yerr=1, uplims=True)

# Aurora's estimate
aurora_hour_offset = 0.17 * 60 #arc second
aurora_hour_freq = 150 # MHz
aurora_hour_scaled = sfaf(aurora_hour_offset, aurora_hour_freq, freq)
legend_lines.append(ax1.plot(freq, aurora_hour_scaled, label="Bulk time offset (1 hour)", color='orange'))
aurora_hour_response = np.exp(-(aurora_hour_scaled / extended_frhm)**2/(2))
ax2.plot(freq, aurora_hour_response, label="Bulk time offset (1 hour)", color='orange')

# Aurora's estimate
aurora_12_offset = 1. * 60 #arc second
aurora_12_freq = 150 # MHz
aurora_12_scaled = sfaf(aurora_12_offset, aurora_12_freq, freq)
legend_lines.append(ax1.plot(freq, aurora_12_scaled, label="Bulk time offset (12 hour)", color='red'))
aurora_12_response = np.exp(-(aurora_12_scaled / extended_frhm)**2/(2))
ax2.plot(freq, aurora_12_response, label="Bulk time offset (12 hour)", color='red')

# Set up plot
# ticks
ax1.set_xscale("log")
ax2.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xticks([80, 95, 115, 140, 170, 205, 250, 300])
ax2.set_xticks([80, 95, 115, 140, 170, 205, 250, 300])
ax1.set_yticks([4, 8, 15, 30, 60, 120, 240])#300, 600, 1200 ])
ax1.minorticks_off()
ax2.minorticks_off()
ax1.get_xaxis().set_major_formatter(ScalarFormatter())
ax2.get_xaxis().set_major_formatter(ScalarFormatter())
ax1.get_yaxis().set_major_formatter(ScalarFormatter())
ax1.set_xlabel("Frequency (MHz)")
ax2.set_xlabel("Frequency (MHz)")
ax1.set_ylabel(r"Anglular Distance ($^{\prime\prime}$)")
ax2.set_ylabel(r"Relative Sensitivty")

#ax1.legend(loc='upper right', ncol=3, bbox_to_anchor=(1.8, 0.88))

fig.legend(legend_lines,     # The line objects
           labels=["P2E HWHM", "",
                   "Residual angular offset (50%)", "Residual angular offset (90%)",
                   "Bulk time offset (1 hour)", "Bulk time offset (12 hour)"],
           loc="upper center",   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box)
           ncol=3,
)
fig.savefig("theoretical_offset.png")
