import find_pulsar_in_obs
from mwa_metadb_utils import getmeta
#obsids = find_pulsar_in_obs.find_obsids_meta_pages(params={'mode':'VOLTAGE_START'})

obs_re = [1118463408, 1114691824, 1113678424, 1113678168, 1107886648, 1102279816, 1102276816, 1102273216, 1102270216, 1102269616, 1101930376, 1101928576, 1101927616, 1101925816, 1101491208, 1100632896, 1100288192, 1099415632, 1099414416, 1098832576, 1098832336, 1098827160, 1098439592, 1097404000, 1096110392, 1095506112]

#Remove obs before  2017-09-30
for obs in obs_re:
    ob_meta = getmeta(params={'obs_id':obs})
    print("Obsid: {} Dur: {:4d}s Project name: {}".format(obs, ob_meta[u'stoptime'] - ob_meta[u'starttime'], ob_meta[u'obsname']))
