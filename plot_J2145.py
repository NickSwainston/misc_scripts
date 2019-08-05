import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import math
from math import sqrt
import numpy as np
import find_pulsar_in_obs as fpio
import mwa_metadb_utils as meta

true_pos_raj  = '21:45:50.46'
true_pos_decj = '-07:50:18.48'
coord = SkyCoord(true_pos_raj, true_pos_decj, unit=(u.hourangle,u.deg))
true_pos_rad = coord.ra.degree 
true_pos_decd = coord.dec.degree


rajs=["22:50:11.54",
"22:15:33.69",
"22:16:38.43",
"22:16:38.51",
"22:19:53.03",
"22:14:28.86",
"22:47:26.01",
"22:52:50.20",
"22:16:06.06",
"22:17:43.33",
"22:17:43.38",
"22:18:15.77",
"22:51:45.38",
"22:15:29.85",
"22:18:48.20",
"22:19:20.59",
"22:15:01.24",
"22:46:55.60",
"22:19:49.53",
"22:45:50.79",
"21:46:53.89",
"21:43:39.27",
"21:44:12.34",
"21:45:18.65",
"21:45:50.46",
"21:45:51.06",
"21:46:23.07",
"22:17:14.11",
"21:45:17.20",
"21:45:18.24",
"21:45:49.01",
"22:16:08.50",
"21:45:50.18",
"21:47:27.90"]

decjs=["-07:11:01.83",
"-02:50:27.39",
"-03:23:55.96",
"-02:50:27.39",
"-02:16:41.21",
"-02:50:27.39",
"+02:43:51.65",
"+03:20:57.88",
"-03:07:13.83",
"-02:50:27.39",
"-02:16:41.21",
"-02:33:36.55",
"+03:20:57.88",
"+06:54:01.86",
"-02:16:41.21",
"-02:33:36.55",
"-03:07:13.83",
"-07:11:01.83",
"+08:57:51.57",
"-07:11:01.83",
"+02:43:51.65",
"+02:43:51.65",
"+02:25:28.57",
"+01:12:53.44",
"-07:50:18.48",
"+01:30:57.56",
"+01:49:01.67",
"-12:03:19.31",
"+02:25:28.57",
"+01:49:01.67",
"+02:43:51.65",
"-12:03:19.31",
"+02:07:11.96",
"+01:49:01.67"]

sns=[7.4,
11.3,
11.5,
9.2,
9.7,
7.0,
8.2,
8.7,
14.3,
14.1,
11.2,
21.8,
7.8,
9.3,
12.9,
8.6,
10.6,
6.7,
8.7,
6.8,
7.0,
9.6,
6.3,
11.4,
166.9,
7.5,
17.0,
23.9,
16.05,
11.5,
15.9,
13.5,
12.6,
6.6]

import matplotlib as mpl
import matplotlib.cm as cm
norm = mpl.colors.Normalize(vmin=min(sns), vmax=max(sns))
cmap = cm.plasma
m = cm.ScalarMappable(norm=norm, cmap=cmap)


coords = SkyCoord(rajs, decjs, unit=(u.hourangle,u.deg))
rads = coords.ra.degree
decds = coords.dec.degree
#fig, axes = plt.subplots(nrows=2)
ax = plt.gca()
#plt.axis('equal')

distance = []
for rad, decd, sn in zip(rads, decds, sns):
    circle = plt.Circle((rad, decd), 0.15, color=m.to_rgba(sn))
    ax.add_patch(circle)
    distance.append(sqrt( (rad - true_pos_rad)**2 + (decd - true_pos_decd)**2))

print("Max distance: {0} degrees".format(max(distance)))
print("J2145-0750 detected {} times".format(len(sns)+1))

#true position
#circle = plt.Circle((true_pos_rad, true_pos_decd), 0.15, color='b')
#ax.add_patch(circle)

#area searched
obs_centre_raj = '22:30:08.79'
obs_centre_decj = '+01:30:57.56'
coord = SkyCoord(obs_centre_raj, obs_centre_decj, unit=(u.hourangle,u.deg))
obs_centre_rad = coord.ra.degree
obs_centre_decd = coord.dec.degree

"""
#calc number of loops
n_beam = 7100
l = 0
n = 1
while n < n_beam:
    n += l * 6
    l += 1
search_radius = 0.9 * 0.3 * float(l)
print("Number of loops: {}".format(l))
circle = plt.Circle((obs_centre_rad, obs_centre_decd), search_radius, color='b', fill=False)
ax.add_patch(circle)
"""


#plot fwhm of tile beam
ob = 1221832280
ob, ra, dec, time, delays,centrefreq, channels = meta.get_common_obs_metadata(ob)
print(delays)
delays = [[12, 12, 12, 12, 8, 8, 8, 8, 4, 4, 4, 4, 0, 0, 0, 0], [12, 12, 12, 12, 8, 8, 8, 8, 4, 4, 4, 4, 0, 0, 0, 0]]
cord = [ob, ra, dec, time, delays, centrefreq, channels]
z=[] ; z_sens =[] ; x=[] ; y=[]

map_dec_range = range(int(dec-15),int(dec+13),1)
map_ra_range = range(int(ra-15),int(ra+37),1)
RA=[] ; Dec=[]
for i in map_dec_range:
    for j in map_ra_range:
        Dec.append(i)
        RA.append(j)

#print(max(Dec), min(RA), Dec.dtype)
time_intervals = 600 # seconds
names_ra_dec = np.column_stack((['source']*len(RA), RA, Dec))
powout = fpio.get_beam_power_over_time(cord, names_ra_dec, dt=time_intervals, degrees = True)
#grab a line of beam power for the pointing declination
#if i == 0:
#    print("len powers list: " + str(powout.shape))
for c in range(len(RA)):
    temppower = 0.
    temppower_sense = 0.
    for t in range(powout.shape[1]):
        power_ra = powout[c,t,0]
        if power_ra > temppower:
            temppower = power_ra
    z.append(temppower)
    x.append(RA[c])
    y.append(Dec[c])
    #x.append(-RA[c]/180.*np.pi +np.pi)
    #y.append(Dec[c]/180.*np.pi)

nx=np.array(x) ; ny=np.array(y); nz=np.array(z)
plt.tricontour(nx, ny, nz, levels=[0.5*max(nz)], alpha = 0.6,
                           colors='blue',
                           linewidths=1.0)

#plot
plt.axis('scaled')
plt.xlabel("Right Acension (degrees)")
plt.ylabel("Declination (degrees)")
gradient = np.linspace(min(sns), max(sns), 256)
gradient = np.vstack((gradient, gradient))
plt.imshow(gradient, origin="right", aspect='auto', cmap=cmap)
plt.colorbar()

#plt.show()
#fig.show()
#plt.colorbar(label="Presto signal to noise")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('J2145-0750_detections.png', dpi=100)
plt.savefig('J2145-0750_detections.eps')
