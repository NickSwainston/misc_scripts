from pulsar_spectra import catalogue
from pulsar_spectra.spectral_fit import find_best_spectral_fit
import psrqpy
import matplotlib.pyplot as plt

query = psrqpy.QueryATNF(loadfromdb=catalogue.ATNF_LOC).pandas

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 3))

#pulsar = 'J0835-4510'
pulsar = 'J0820-1350'

antf_dict = catalogue.all_flux_from_atnf(query=query)
# Resort into cat_list format
antf_cat = { pulsar : [ [], [], [], [] ] }
for ref in antf_dict[pulsar].keys():
    # Update list
    for freq, flux, flux_err in zip(antf_dict[pulsar][ref]['Frequency MHz'],
                                    antf_dict[pulsar][ref]['Flux Density mJy'],
                                    antf_dict[pulsar][ref]['Flux Density error mJy']):
        antf_cat[pulsar][0] += [freq]
        antf_cat[pulsar][1] += [flux]
        antf_cat[pulsar][2] += [flux_err]
        antf_cat[pulsar][3] += [ref]
freqs, fluxs, flux_errs, refs = antf_cat[pulsar]
model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, freqs, fluxs, flux_errs, refs, plot_best=True, alternate_style=True, axis=axs[0])
#axs[ax_i//cols, ax_i%cols].set_title('PSR '+pulsar)


cat_dict = catalogue.collect_catalogue_fluxes(query=query)
# reorganise so it is in the same order
freqs, fluxs, flux_errs, refs = cat_dict[pulsar]
fo = []
flo = []
fleo = []
ro = []
for ref in antf_dict[pulsar].keys():
    for fi, fli, flei, ri in zip(freqs, fluxs, flux_errs, refs):
        if ref == ri:
            fo.append(fi)
            flo.append(fli)
            fleo.append(flei)
            ro.append(ri)
# put the rest in
for fi, fli, flei, ri in zip(freqs, fluxs, flux_errs, refs):
    if ri not in antf_dict[pulsar].keys():
        fo.append(fi)
        flo.append(fli)
        fleo.append(flei)
        ro.append(ri)
model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, fo, flo, fleo, ro, plot_best=True, alternate_style=True, axis=axs[1])
#axs[ax_i//cols, ax_i%cols].set_title('PSR '+pulsar)

plt.tight_layout(pad=2.5)
plt.savefig(f"antf_comparison_{pulsar}.png", bbox_inches='tight', dpi=300)
