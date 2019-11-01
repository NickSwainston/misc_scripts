import json
import urllib.request

#download the json and turn into a list of each SN which each has a dictionary of all it's values
sn_data = json.load(urllib.request.urlopen('https://raw.githubusercontent.com/astrocatalogs/supernovae/master/output/catalog.json'))

# loop over the SNs and only keep the Ias
Ia_sns = []
for sn in sn_data:
    #get all the keys which are basically the variables this SN has
    snkeys = sn.keys()
    if 'claimedtype' in snkeys:
        if sn['claimedtype'][0]['value'] == 'Ia':
            #if Ia record it in a new list
            Ia_sns.append(sn)

for sn in Ia_sns:
    snkeys = sn.keys()
    if 'ra' in snkeys:
        ra = sn['ra'][0]['value']
    else:
        ra = None
    if 'dec' in snkeys:
        dec = sn['dec'][0]['value']
    else:
        dec = None
    if 'name' in snkeys:
        name = sn['name']
    else:
        name = None
    print('{} {} {}'.format(name, ra, dec))
