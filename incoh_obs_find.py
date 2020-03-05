import find_pulsar_in_obs
obsids = find_pulsar_in_obs.find_obsids_meta_pages(params={'mode':'VOLTAGE_START', 'contigfreq':1, 'cenchan':145, 'dataquality':1})

#Remove obs before  2017-09-30
bobsids = []
for obs in obsids:
    if int(obs) < 1190762640:
        bobsids.append(obs)

print('Number of obs: {}'.format(len(bobsids)))
print(bobsids)


