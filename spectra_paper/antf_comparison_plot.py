from pulsar_spectra import catalogue
from pulsar_spectra.spectral_fit import find_best_spectral_fit
import matplotlib.pyplot as plt

#pulsars = ['J0820-1350', 'J0835-4510', 'J1820-0427']
#pulsars = ['J0835-4510', 'J0953+0755']
pulsars = [
    # BPL
    # 'J0835-4510',
    'J0543+2329',
    # LFTO
    # 'J1751-4657',
    # 'J1509+5531',
    'J1932+1059',
    # HFTO
    # 'J0953+0755',
    # 'J1012+5307',
    # 'J1740+1311',
    'J1835-0643',
    # DTOS
    # 'J1607-0032'
    # 'J1932+1059',
    'J0837+0610',
]
np = len(pulsars)

fig, axs = plt.subplots(nrows=np, ncols=2, figsize=(10, 4*np))

antf_dict = catalogue.all_flux_from_atnf()
# Resort into cat_list format
for pi, pulsar in enumerate(pulsars):
    antf_refs = []
    antf_cat = { pulsar : [ [], [], [], [], [] ] }
    for ref in antf_dict[pulsar].keys():
        if ref in ("Mignani_2017_ATNF", "Taylor_1993_ATNF"):
            continue
        antf_refs.append(ref[:-5])
        # Update list
        for freq, flux, flux_err in zip(antf_dict[pulsar][ref]['Frequency MHz'],
                                        antf_dict[pulsar][ref]['Flux Density mJy'],
                                        antf_dict[pulsar][ref]['Flux Density error mJy']):
            antf_cat[pulsar][0] += [freq]
            antf_cat[pulsar][1] += [None]
            antf_cat[pulsar][2] += [flux]
            antf_cat[pulsar][3] += [flux_err]
            antf_cat[pulsar][4] += [ref]
    freqs, bands, fluxs, flux_errs, refs = antf_cat[pulsar]
    model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, freqs, bands, fluxs, flux_errs, refs, plot_best=True, alternate_style=True, axis=axs[pi][0])
    axs[pi][0].set_title(f'PSR {pulsar} ATNF')


    cat_dict = catalogue.collect_catalogue_fluxes()
    # reorganise so it is in the same order
    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    fo = []
    bo = []
    flo = []
    fleo = []
    ro = []
    for ref in antf_refs:
        for fi, bi, fli, flei, ri in zip(freqs, bands, fluxs, flux_errs, refs):
            if ref == ri:
                fo.append(fi)
                bo.append(bi)
                flo.append(fli)
                fleo.append(flei)
                ro.append(ri)
    # put the rest in
    print(refs)
    print(antf_refs)
    for fi, bi, fli, flei, ri in zip(freqs, bands, fluxs, flux_errs, refs):
        if ri not in antf_refs:
            fo.append(fi)
            bo.append(bi)
            flo.append(fli)
            fleo.append(flei)
            ro.append(ri)

    for freq, band, flux, flux_err, ref in zip(fo, bo, flo, fleo, ro):
        print(f"{str(freq):8s}{str(band):8s}{float(flux):8.2f}{float(flux_err):8.2f} {str(ref):20s}")
    model, m, fit_info, p_best, p_category = find_best_spectral_fit(pulsar, fo, bo, flo, fleo, ro, plot_best=True, alternate_style=True, axis=axs[pi][1])
    axs[pi][1].set_title(f'PSR {pulsar} pulsar_spectra')
    print(model)

    #make same y axis the same
    print(fluxs)
    ymin = min(flo + fluxs) * 0.1
    ymax = max(flo + fluxs) ** 1.1
    axs[pi][0].set_ylim(ymin, ymax)
    axs[pi][1].set_ylim(ymin, ymax)
    #make same x axis the same
    print(fo)
    xmin = min(fo) * 0.7
    xmax = max(fo) * 1.5
    axs[pi][0].set_xlim(xmin, xmax)
    axs[pi][1].set_xlim(xmin, xmax)
#axs[ax_i//cols, ax_i%cols].set_title('PSR '+pulsar)

plt.tight_layout(pad=2.5)
pulsar_str = "_".join(pulsars)
plt.savefig(f"antf_comparison_{pulsar_str}.png", bbox_inches='tight', dpi=300)
