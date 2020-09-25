import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

use_sns = True
if use_sns:
    fig, ax = plt.subplots(1, 3, figsize=(10, 5))
else:
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))


#FHRM radius
p2cr = 18.56 / 60. /2.
p1r  =  2.53 / 60. /2.
p2er = p1r / 2

#1255444104 20 mins data, pdmp SN
p2c_grid = [[8.75, "00:35:25.20_-10:35:14.03"],
            [11.59, "00:35:44.07_-10:43:56.49"],
            [11.09, "00:35:44.09_-10:26:31.56"],
            [21.67, "00:36:02.97_-10:35:14.03"],
            #[22.61, "00:36:09.74_-10:32:11.18"],
            [16.67, "00:36:21.84_-10:26:31.56"],
            [13.07, "00:36:21.86_-10:43:56.49"],
            [14.29, "00:36:40.73_-10:35:14.03"]]

p2c_ras = []
p2c_decs = []
p2c_sns = []
for g in p2c_grid:
    ra, dec = g[1].split('_')
    p2c_ras.append(ra)
    p2c_decs.append(dec)    
    p2c_sns.append(g[0])

#p2c grid
coord = SkyCoord(p2c_ras, p2c_decs, unit=(u.hourangle,u.deg))
p2c_decs_deg = coord.dec.degree
p2c_ras_deg = coord.ra.degree
for ra, dec, sn in zip(p2c_ras_deg, p2c_decs_deg, p2c_sns):
    if ra in [min(p2c_ras_deg), max(p2c_ras_deg)]:
        linewidth = 4
    else:
        linewidth = 1
    if use_sns:
        ax[0].text(ra, dec, str(sn), color='g', fontsize=8, ha='center', va='center')
        circle = plt.Circle((ra, dec), p2cr, color='g', fill=False, linewidth=linewidth, alpha=sn/max(p2c_sns))
    else:
        circle = plt.Circle((ra, dec), p2cr, color='g', fill=False, linewidth=linewidth)
    ax[0].add_artist(circle)


#p2c posititon estimate
# old random estimate "00:36:14.25", "-10:35:14.03"
coord = SkyCoord("00:36:09.74", "-10:32:11.18", unit=(u.hourangle,u.deg))
dec_est = coord.dec.degree
ra_est = coord.ra.degree
circle = plt.Circle((ra_est, dec_est), p2cr/2, color='g', alpha=0.5)
ax[0].add_artist(circle)

ax[0].set_ylim(p2c_decs_deg[0]-2*p2cr, p2c_decs_deg[0]+2*p2cr)
ax[0].set_xlim(min(p2c_ras_deg) - p2cr, max(p2c_ras_deg) + p2cr)
ax[0].set_xlabel("Right Acension (degrees)")
ax[0].set_ylabel("Declination (degrees)")
#plt.gca().set_aspect('equal', adjustable='box')
ax[0].set_aspect('equal', adjustable='box')


#zoom in on the second subplot
circle = plt.Circle((ra_est, dec_est), p2cr/2, color='g', fill=False)
ax[1].add_artist(circle)


p1_grid = [["5.32", "00:35:58.79_-10:35:14.03"],
           ["5.21", "00:36:01.36_-10:34:02.73"],
           ["5.13", "00:36:01.36_-10:36:25.32"],
           ["5.54", "00:36:03.94_-10:32:51.42"],
           ["5.39", "00:36:03.94_-10:35:14.03"],
           ["5.19", "00:36:03.94_-10:37:36.60"],
           ["5.24", "00:36:06.51_-10:38:47.87"],
           ["5.20", "00:36:06.52_-10:31:40.10"],
           ["5.46", "00:36:06.52_-10:34:02.73"],
           ["5.50", "00:36:06.52_-10:36:25.32"],
           ["5.37", "00:36:09.09_-10:32:51.42"],
           ["5.98", "00:36:09.09_-10:35:14.03"],
           ["5.50", "00:36:09.09_-10:37:36.60"],
           ["5.38", "00:36:11.67_-10:31:40.10"],
           ["5.69", "00:36:11.67_-10:34:02.73"],
           ["5.54", "00:36:11.67_-10:36:25.32"],
           ["4.92", "00:36:11.67_-10:38:47.87"],
           ["5.52", "00:36:14.25_-10:32:51.42"],
           ["5.44", "00:36:14.25_-10:35:14.03"],
           ["5.38", "00:36:14.25_-10:37:36.60"],
           ["6.16", "00:36:16.82_-10:31:40.10"],
           ["5.46", "00:36:16.82_-10:34:02.73"],
           ["5.25", "00:36:16.82_-10:36:25.32"],
           ["5.18", "00:36:16.82_-10:38:47.87"],
           #["5.74", "00:36:17.70_-10:34:04.20"],
           ["5.14", "00:36:19.40_-10:32:51.42"],
           ["5.42", "00:36:19.40_-10:35:14.03"],
           ["5.20", "00:36:19.40_-10:37:36.60"],
           ["4.92", "00:36:21.97_-10:31:40.10"],
           ["5.37", "00:36:21.97_-10:34:02.73"],
           ["5.54", "00:36:21.98_-10:36:25.32"],
           ["5.52", "00:36:21.98_-10:38:47.87"],
           ["5.08", "00:36:24.55_-10:32:51.42"],
           ["5.36", "00:36:24.55_-10:35:14.03"],
           ["5.02", "00:36:24.55_-10:37:36.60"],
           ["6.30", "00:36:27.13_-10:34:02.73"],
           ["5.15", "00:36:27.13_-10:36:25.32"],
           ["5.37", "00:36:29.70_-10:35:14.03"]]

for pos in p1_grid:
    raj, decj = pos[1].split("_")
    coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
    dec = coord.dec.degree
    ra = coord.ra.degree
    circle = plt.Circle((ra, dec), p1r, color='r', fill=False, linewidth=0.5)
    ax[1].add_artist(circle)

p2e_grid_final = [[17.67, "00:36:14.70_-10:33:37.18"],
                  [8.59, "00:36:12.13_-10:33:37.18"],
                  [5.21, "00:36:13.41_-10:33:01.70"],
                  [21.51, "00:36:15.98_-10:33:01.70"],
                  [21.60, "00:36:17.26_-10:33:37.18"],
                  [13.65, "00:36:15.98_-10:34:12.65"],
                  [8.95, "00:36:13.41_-10:34:12.65"]]

p2e_sns = []
for sn in p2e_grid_final:
    p2e_sns.append(sn[0])

for pos in p2e_grid_final:
    raj, decj = pos[1].split("_")
    coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
    dec = coord.dec.degree
    ra = coord.ra.degree
    sn = pos[0]
    if use_sns:
        ax[1].text(ra, dec, str(sn), color='b', fontsize=8, ha='center', va='center')
        ax[2].text(ra, dec, str(sn), color='b', fontsize=8, ha='center', va='center')
        circle = plt.Circle((ra, dec), p2er, color='b', fill=False, linewidth=0.5, alpha=sn/max(p2e_sns))
        ax[2].add_artist(circle)
        circle = plt.Circle((ra, dec), p2er, color='b', fill=False, linewidth=0.5, alpha=sn/max(p2e_sns))
    else:
        circle = plt.Circle((ra, dec), p2er, color='b', fill=False, linewidth=0.5)
    ax[1].add_artist(circle)

ax[1].set_ylim(dec_est-p2cr/2, dec_est+p2cr/2)
ax[1].set_xlim(ra_est-p2cr/2, ra_est+p2cr/2)
ax[1].set_xlabel("Right Acension (degrees)")
ax[1].yaxis.tick_right()
ax[1].set_aspect('equal', adjustable='box')

if use_sns:
    coord = SkyCoord("00:36:14.70", "-10:33:37.18", unit=(u.hourangle,u.deg))
    dec = coord.dec.degree
    ra = coord.ra.degree
    ax[2].set_ylim(dec-2*p2er, dec+2*p2er)
    ax[2].set_xlim(ra-2*p2er, ra+2*p2er)
    ax[2].set_xlabel("Right Acension (degrees)")
    ax[2].yaxis.tick_right()
    ax[2].set_aspect('equal', adjustable='box')

coord = SkyCoord("00:36:14.61", "-10:33:42.81", unit=(u.hourangle,u.deg))
dec = coord.dec.degree
ra = coord.ra.degree
circle = plt.Circle((ra, dec), p2er/2, color='b', alpha=0.5)
ax[1].add_artist(circle)

#mark zoomed in subplot
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax[0], ax[1], loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)
if use_sns:
    mark_inset(ax[1], ax[2], loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)

plt.savefig("first_mwa_discovery_localisation.png", dpi=1000)