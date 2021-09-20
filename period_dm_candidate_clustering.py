import csv
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky

with open('data/1226062160_positive_detections_list.txt') as csvfile:
#with open('data/1226062160_all_detections_list.txt') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ')
    data = []
    ra = []
    dec = []
    period = []
    dm = []
    sn = []
    for row in spamreader:
        #RA Dec Period (ms) DM
        data.append(row)
        ra.append(row[0])
        dec.append(row[1])
        period.append(float(row[2]))
        dm.append(float(row[3]))
        sn.append(float(row[4]))

#RA Dec Period(ms) DM
data = np.array(data)
period = np.array(period)
dm = np.array(dm)
sn = np.array(sn)

# Maximum differences
mdist = 18.56 # arcminutes
mperiod = 1 # ms
mdm = 5

files_to_follow_up = []

print("Making coords")
coords = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
print("Searching")
# Check if they're within a beam width
idx1, idx2, sep2d, _ = search_around_sky(coords, coords, mdist*u.arcminute)

period = []
dm = []
sn = []
for i in range(len(sep2d)):
    cand_1 = data[idx1[i]]
    cand_2 = data[idx2[i]]
    if ( idx1[i] != idx2[i]) and \
        ( abs(float(cand_1[2]) - float(cand_2[2])) < mperiod ) and \
        ( abs(float(cand_1[3]) - float(cand_2[3])) < mdm ) and \
        not ( 1642. < float(cand_1[2]) < 1644. ):
        print("")
        print(cand_1)
        print(cand_2)
        period.append(float(cand_1[2]))
        dm.append(float(cand_1[3]))
        sn.append(float(cand_1[4]))
        files_to_follow_up.append(cand_1[5][:-8]+"png")
        files_to_follow_up.append(cand_2[5][:-8]+"png")
files_to_follow_up = pd.unique(files_to_follow_up).tolist()

#sort by SN
period = np.array(period)
dm = np.array(dm)
sn = np.array(sn)
idx = sn.argsort()
dm, period, sn = dm[idx], period[idx], sn[idx]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(dm, period, s=5, c=sn, norm=LogNorm())
#plt.scatter(sn, dm)
ax.set_yscale('log')
plt.show()

with open('clustered_candidates.txt', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for file_name in files_to_follow_up:
            spamwriter.writerow([file_name])
