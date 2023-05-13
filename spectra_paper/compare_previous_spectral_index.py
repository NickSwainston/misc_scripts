from pulsar_spectra.catalogue import collect_catalogue_fluxes, CAT_DIR
import yaml
#cat_dict = collect_catalogue_fluxes()

previous_publications = [
    "Spiewak_2022",
    "Jankowski_2018",
    #"Han_2017",
    "Bilous_2016",
    #"Maron_2000",
    "Malofeev_2000",
    "Toscano_1998",
    #"Lorimer_1995",
]

for pub in previous_publications:
    with open(f"{CAT_DIR}/{pub}.yaml", "r") as stream:
        cat_dict = yaml.safe_load(stream)
    pulsars = cat_dict.keys()
    print(f"\cite{{{pub.replace('_', ''):20s}}} & {len(pulsars)} \\\\")

