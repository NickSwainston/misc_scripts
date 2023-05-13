import numpy as np
import matplotlib.pyplot as plt

from pulsar_spectra.spectral_fit import find_best_spectral_fit
from pulsar_spectra.catalogue import collect_catalogue_fluxes



cat_dict = collect_catalogue_fluxes()
pulsars = [
    # Likely true
    "J0024-7204C",
    "J2124-3358",

    # Not enough points as it curves
    "J0024-7204J",
    "J0621+1002",

    # Curves lower than reasonable
    "J0030+0451",
    "J1300+1240",

    # ignores some low frequency points
    "J1600-3053",
    "J1623-2631",
]
fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(10, 4*3))

alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']

# Resort into cat_list format
for pi, pulsar in enumerate(pulsars):
    yi = pi % 2
    xi = pi // 2
    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, freqs, bands, fluxs, flux_errs, refs, plot_best=True, alternate_style=True, axis=axs[xi][yi])
    axs[xi][yi].set_title(f"({alphabet[pi]}) {pulsar}")
    print(model)

#axs[ax_i//cols, ax_i%cols].set_title('PSR '+pulsar)

plt.tight_layout(pad=2.5)
plt.savefig(f"msp_turn_over.png", bbox_inches='tight', dpi=300)
