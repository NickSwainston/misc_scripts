import numpy as np
import matplotlib.pyplot as plt

from pulsar_spectra.spectral_fit import find_best_spectral_fit
from pulsar_spectra.catalogue import collect_catalogue_fluxes



cat_dict = collect_catalogue_fluxes()
pulsars = [
    ('J0415+6954', "Simple power law"),
    ('J0437-4715', "Broken power law"),
    ('J0809-4753', "Low-frequency turn over power law"),
    ('J1703-1846', "High-frequency cut off power law"),
    ('J0946+0951', "Log parabolic spectrum"),
    ('J1932+1059', "Double turn over spectrum"),
]
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10, 4*3))

# Resort into cat_list format
for pi, pulsar_pair in enumerate(pulsars):
    pulsar, exp_model_name = pulsar_pair
    yi = pi % 2
    xi = pi // 2
    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    if pulsar == 'J1932+1059':
        # Add an extra point to encourage a high frequency cut off
        freqs.append(25000)
        bands.append(1000)
        fluxs.append(0.2)
        flux_errs.append(0.02)
        refs.append("Bartel_1978")
    model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, freqs, bands, fluxs, flux_errs, refs, plot_best=True, alternate_style=True, axis=axs[xi][yi])
    axs[xi][yi].set_title(exp_model_name)
    print(model)

#axs[ax_i//cols, ax_i%cols].set_title('PSR '+pulsar)

plt.tight_layout(pad=2.5)
plt.savefig(f"model_examples.png", bbox_inches='tight', dpi=300)
