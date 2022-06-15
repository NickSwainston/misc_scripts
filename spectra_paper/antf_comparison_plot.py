from pulsar_spectra import catalogue
from pulsar_spectra.spectral_fit import find_best_spectral_fit
import psrqpy
import matplotlib.pyplot as plt

query = psrqpy.QueryATNF(loadfromdb=catalogue.ATNF_LOC).pandas
#pulsars = ['J0820-1350', 'J0835-4510', 'J1820-0427']
pulsars = ['J0835-4510', 'J1141-6545', 'J1751-4657', 'J0953+0755']
np = len(pulsars)

fig, axs = plt.subplots(nrows=np, ncols=2, figsize=(9, 3.5*np))

antf_dict = catalogue.all_flux_from_atnf(query=query)
# Resort into cat_list format
for pi, pulsar in enumerate(pulsars):
    antf_cat = { pulsar : [ [], [], [], [] ] }
    for ref in antf_dict[pulsar].keys():
        if ref == "Mignani_2017_ATNF":
            continue
        # Update list
        for freq, flux, flux_err in zip(antf_dict[pulsar][ref]['Frequency MHz'],
                                        antf_dict[pulsar][ref]['Flux Density mJy'],
                                        antf_dict[pulsar][ref]['Flux Density error mJy']):
            antf_cat[pulsar][0] += [freq]
            antf_cat[pulsar][1] += [flux]
            antf_cat[pulsar][2] += [flux_err]
            antf_cat[pulsar][3] += [ref]
    freqs, fluxs, flux_errs, refs = antf_cat[pulsar]
    model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, freqs, fluxs, flux_errs, refs, plot_best=True, alternate_style=True, axis=axs[pi][0])
    axs[pi][0].set_title(f'PSR {pulsar} ATNF')


    cat_dict = catalogue.collect_catalogue_fluxes(query=query)
    # reorganise so it is in the same order
    freqs, fluxs, flux_errs, refs = cat_dict[pulsar]
    fo = []
    flo = []
    fleo = []
    ro = []
    # for ref in antf_dict[pulsar].keys():
    #     for fi, fli, flei, ri in zip(freqs, fluxs, flux_errs, refs):
    #         if ref[:-5] == ri:
    #             fo.append(fi)
    #             flo.append(fli)
    #             fleo.append(flei)
    #             ro.append(ri)
    # put the rest in
    for fi, fli, flei, ri in zip(freqs, fluxs, flux_errs, refs):
        if ri not in antf_dict[pulsar].keys():
            fo.append(fi)
            flo.append(fli)
            fleo.append(flei)
            ro.append(ri)

    for freq, flux, flux_err, ref in zip(fo, flo, fleo, ro):
        print(f"{str(freq):8s}{float(flux):8.2f}{float(flux_err):8.2f} {str(ref):20s}")
    model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, fo, flo, fleo, ro, plot_best=True, alternate_style=True, axis=axs[pi][1])
    axs[pi][1].set_title(f'PSR {pulsar} pulsar_spectra')

    #make same y axis the same
    print(fluxs)
    ymin = min(flo + fluxs) * 0.1
    ymax = max(flo + fluxs) ** 1.1
    axs[pi][0].set_ylim(ymin, ymax)
    axs[pi][1].set_ylim(ymin, ymax)
#axs[ax_i//cols, ax_i%cols].set_title('PSR '+pulsar)

plt.tight_layout(pad=2.5)
pulsar_str = "_".join(pulsars)
plt.savefig(f"antf_comparison_{pulsar_str}.png", bbox_inches='tight', dpi=300)
