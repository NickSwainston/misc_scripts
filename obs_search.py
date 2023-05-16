from vcstools.metadb_utils import getmeta, find_obsids_meta_pages

params = {'mode':'VOLTAGE_START', 'anychan':63}
obsids= find_obsids_meta_pages(params)

for obsid in obsids:
    if get_obs_array_phase(obsid) == 'P2E':
        obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels= get_common_obs_metadata(obsid)
        print(obisd, dura, centrefreq)
