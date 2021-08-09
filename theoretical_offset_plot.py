
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

fig, ax = plt.subplots()
#ax.plot(freq, compact_frhm, label="P2C FRHM")
ax.plot(freq, extended_frhm, label="P2E FRHM")

sfaf = np.vectorize(scale_for_all_freqs)
# Chris's estimate
chris_offset = 0.35 * 60 #arc second
chris_freq = 200 # MHz
chris_scaled = sfaf(chris_offset, chris_freq, freq)
ax.plot(freq, chris_scaled, label="Chris estimate")

#plt.errorbar(freq, chris_scaled, yerr=1, uplims=True)

# Aurora's estimate
aurora_offset = 1. * 60 #arc second
aurora_freq = 150 # MHz
aurora_scaled = sfaf(aurora_offset, aurora_freq, freq)
ax.plot(freq, aurora_scaled, label="Aurora estimate")

# Set up plot
# ticks
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([80, 95, 115, 140, 170, 205, 250, 300])
ax.set_yticks([15, 30, 60, 120, 240])#300, 600, 1200 ])
plt.minorticks_off()
ax.get_xaxis().set_major_formatter(ScalarFormatter())
ax.get_yaxis().set_major_formatter(ScalarFormatter())
ax.set_xlabel("Frequency (MHz)")
ax.set_ylabel(r"Angle ($^{\prime\prime}$)")

fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.88))
fig.savefig("theoretical_offset.png")
