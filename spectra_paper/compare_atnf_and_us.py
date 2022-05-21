import psrqpy
import glob
import yaml

from pulsar_spectra.catalogue import all_flux_from_atnf,  collect_catalogue_fluxes



cat_dict = collect_catalogue_fluxes()

atnf = []
both = []
us = []
none = []
for jname in cat_dict.keys():
    atnf_bool = False
    us_bool = False
    for ref in cat_dict[jname][3]:
        if "ATNF" in ref:
            atnf_bool = True
        else:
            us_bool = True

    if atnf_bool and us_bool:
        both.append(jname)
    elif atnf_bool:
        atnf.append(jname)
    elif us_bool:
        us.append(jname)
    else:
        none.append(jname)

print(f"atnf: {len(atnf)}")
print(f"both: {len(both)}")
print(f"us  : {len(us)}")
print(f"none: {len(none)}")