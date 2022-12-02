from pulsar_spectra.catalogue import collect_catalogue_fluxes



cat_dict = collect_catalogue_fluxes()

atnf = []
both = []
us = []
none = []
for jname in cat_dict.keys():
    atnf_bool = False
    us_bool = False
    for ref in cat_dict[jname][-1]:
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


ref_count = {}
for pulsar in cat_dict.keys():
    refs = cat_dict[pulsar][-1]
    for ref in refs:
        if 'ATNF' in ref:
            if ref in ref_count.keys():
                ref_count[ref] += 1
            else:
                ref_count[ref] = 1
count_list = []
for ref in ref_count.keys():
    count_list.append([ref, ref_count[ref]])
count_list.sort(key = lambda x: x[1], reverse=True)
for ref, num in count_list[:20]:
    print(f"{ref:25s} {num}")
