import matplotlib.pyplot as plt
from pulsar_spectra.spectral_fit import find_best_spectral_fit
from pulsar_spectra.catalogue import collect_catalogue_fluxes

# Pulsar, flux, flux_err
pulsars = [
    'J1012+5307',
    'J1740+1311',
]
cols = 1
rows = 2
fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=(5*cols, 4*rows))

cat_list = collect_catalogue_fluxes()
for pi, pulsar in enumerate(pulsars):
    freqs, bands, fluxs, flux_errs, refs = cat_list[pulsar]
    model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, freqs, bands, fluxs, flux_errs, refs, plot_best=True, alternate_style=True, axis=axs[pi])
    axs[pi].set_title('PSR '+pulsar)

plt.tight_layout(pad=2.5)
plt.savefig("hfco_examples.png", bbox_inches='tight', dpi=300)