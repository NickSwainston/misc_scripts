
from pulsar_spectra import catalogue


antf_dict = catalogue.all_flux_from_atnf()
# Resort into cat_list format
n_pulsar = 0
for pulsar in antf_dict.keys():
    freqs = []
    for ref in antf_dict[pulsar].keys():
        # Update list
        for freq in antf_dict[pulsar][ref]['Frequency MHz']:
            freqs.append(freq)
    if len(freqs) == 0:
        continue
    if min(freqs) < 300:
        n_pulsar += 1

print(f"ANTF: {n_pulsar} below 300 MHz")


cat_dict = catalogue.collect_catalogue_fluxes()
# Resort into cat_list format
n_pulsar = 0
for pulsar in cat_dict.keys():
    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    if len(freqs) == 0:
        continue
    if min(freqs) < 300:
        n_pulsar += 1

print(f"PS:   {n_pulsar} below 300 MHz")

# Resort into cat_list format
n_pulsar = 0
for pulsar in cat_dict.keys():
    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    if len(freqs) < 4:
        continue
    if min(freqs) < 300:
        n_pulsar += 1


print(f"PS:   {n_pulsar} below 300 MHz")