import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable

import mimic_alpha as ma

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits

import numpy as np

use_sns = True

#fig = plt.figure(constrained_layout=True, figsize=(15, 6))
fig = plt.figure(figsize=(15, 6))
"""
widths = [3, 3, 3, 1]
heights = [1, 3]
#widths = [1, 1, 1, 0.3]
#heights = [0.3, 1]
spec = fig.add_gridspec(ncols=4, nrows=2,
                        width_ratios=widths, height_ratios=heights,
                        hspace=0.)
#gs_kw = dict(width_ratios=widths, height_ratios=heights)
#fig, ax = plt.subplots(figsize=(10, 4),
#                       ncols=4, nrows=2,
#                       constrained_layout=True, gridspec_kw=gs_kw)
fig.add_gridspec(4, 10)

#gs = fig.add_gridspec(4, 10)

p2c_ax = fig.add_subplot(spec[1, 0])
p1_ax  = fig.add_subplot(spec[1, 1])
p2e_ax = fig.add_subplot(spec[1, 2])
x_ax = fig.add_subplot(spec[0, 2])#, sharex = p2e_ax)
x_ax.set_xticks([])
x_ax.set_yticks([])
y_ax = fig.add_subplot(spec[1, 3])#, sharey = p2e_ax)
y_ax.set_xticks([])
y_ax.set_yticks([])
p2e_ax.get_shared_x_axes().join(p2e_ax, x_ax)
p2e_ax.get_shared_y_axes().join(p2e_ax, y_ax)

#blank1 = fig.add_subplot(spec[0, 0])
#blank2 = fig.add_subplot(spec[0, 1])
#blank3 = fig.add_subplot(spec[0, 3])
"""
"""
gs1 = fig.add_gridspec(nrows=4, ncols=10)
p2c_ax = fig.add_subplot(gs1[1:-1,:3])
p1_ax  = fig.add_subplot(gs1[1:-1,3:6])

gs2 = fig.add_gridspec(nrows=4, ncols=10,
                       hspace=0., wspace=0.)
#p2e_ax = fig.add_subplot(gs2[1:-1,:3])
p2e_ax = fig.add_subplot(gs2[1:-1,6:9])
#x_ax = fig.add_subplot(gs2[0,:3])#, sharex = p2e_ax)
x_ax = fig.add_subplot(gs2[0,6:9])#, sharex = p2e_ax)
#y_ax = fig.add_subplot(gs2[1:-1,3])#, sharey = p2e_ax)
y_ax = fig.add_subplot(gs2[1:-1,9])#, sharey = p2e_ax)
"""
#grid = ImageGrid(fig, 111,  # similar to subplot(111)
#                 nrows_ncols=(1, 1),
#                 axes_pad=0.2,
#                 add_all=True,
#                 label_mode="L",
#                 )
#p2c_ax = grid[0]
fig, (p2c_ax, p1_ax, p2e_ax) = plt.subplots(1,3, figsize=(18, 6))
divider = make_axes_locatable(p2e_ax)

#p1_ax = divider.append_axes("left", size=3., pad=0.2, sharex=grid[0])
#p2e_ax = divider.append_axes("left", size=3., pad=3.2, sharex=grid[0])
#p1_ax  = grid[1]
#p2e_ax = grid[2]

fig.savefig('test_grid_first_3.png')


#divider = make_axes_locatable(p2e_ax)
x_ax = divider.append_axes("top", size=0.8, pad=0., sharex=p2e_ax)
y_ax = divider.append_axes("right", size=0.8, pad=0., sharey=p2e_ax)

y_ax.set_xticks([])
#y_ax.set_yticks([])
#x_ax.set_xticks([])
x_ax.set_yticks([])
plt.setp(x_ax.get_xticklabels(), visible=False)
plt.setp(y_ax.get_yticklabels(), visible=False)

fig.savefig('test_grid.png')

#FWHM radius
p2cr = 18.56 / 60. /2.
p1r  =  2.53 / 60. /2.
p2er = p1r / 2

# Phase 2 compact plot --------------------------------------------------------------------
#1255444104 20 mins data, pdmp SN
p2c_grid = [[8.7, "00:35:25.20_-10:35:14.03"],
            [11.5, "00:35:44.07_-10:43:56.49"],
            [11.0, "00:35:44.09_-10:26:31.56"],
            [21.6, "00:36:02.97_-10:35:14.03"],
            #[22.6, "00:36:09.74_-10:32:11.18"],
            [16.6, "00:36:21.84_-10:26:31.56"],
            [13.0, "00:36:21.86_-10:43:56.49"],
            [14.2, "00:36:40.73_-10:35:14.03"]]

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
    if use_sns:
        p2c_ax.text(ra, dec, str(sn), color='black', fontsize=8, ha='center', va='center')
        if ra == min(p2c_ras_deg) or ra == max(p2c_ras_deg):
            circle = plt.Circle((ra, dec), p2cr, color='g', fill=False, linewidth=3)#)#, alpha=sn/max(p2c_sns))
        else:
            circle = plt.Circle((ra, dec), p2cr, color='g', fill=False, linewidth=1)
    else:
        circle = plt.Circle((ra, dec), p2cr, color='g', fill=False, linewidth=1)
    p2c_ax.add_artist(circle)


#p2c posititon estimate
# old random estimate "00:36:14.25", "-10:35:14.03"
coord = SkyCoord("00:36:09.74", "-10:32:11.18", unit=(u.hourangle,u.deg))
dec_est = coord.dec.degree
ra_est = coord.ra.degree
#circle = plt.Circle((ra_est, dec_est), p2cr/2, color='g', alpha=0.5)
circle = plt.Circle((ra_est, dec_est), p2cr/2, color=ma.colorAlpha_to_rgb('g', 0.5)[0], zorder=0.5)
p2c_ax.add_artist(circle)

p2c_ax.set_ylim(p2c_decs_deg[0]-2*p2cr, p2c_decs_deg[0]+2*p2cr)
p2c_ax.set_xlim(min(p2c_ras_deg) - p2cr, max(p2c_ras_deg) + p2cr)
#plt.gca().set_aspect('equal', adjustable='box')
p2c_ax.set_aspect('equal', adjustable='box')


# Phase 1 plot --------------------------------------------------------------------
#zoom in on the second subplot
circle = plt.Circle((ra_est, dec_est), p2cr/2, color='g', fill=False)
p1_ax.add_artist(circle)


p1_grid = [[7.14, "00:36:09.74_-10:32:11.18"],
[5.36, "00:36:04.58_-10:32:11.18"],
[6.01, "00:36:07.16_-10:30:59.87"],
[5.26, "00:36:12.31_-10:30:59.87"],
[6.0, "00:36:14.89_-10:32:11.18"],
[11.3, "00:36:12.31_-10:33:22.48"],
[5.44, "00:36:07.16_-10:33:22.48"],
[5.29, "00:35:59.43_-10:32:11.18"],
[5.92, "00:36:02.01_-10:30:59.87"],
[5.54, "00:36:04.58_-10:29:48.54"],
[6.02, "00:36:09.74_-10:29:48.54"],
[5.59, "00:36:14.89_-10:29:48.54"],
[5.13, "00:36:17.46_-10:30:59.87"],
[5.56, "00:36:20.04_-10:32:11.18"],
[17.3, "00:36:17.46_-10:33:22.48"],
[19.1, "00:36:14.89_-10:34:33.78"],
[5.31, "00:36:09.74_-10:34:33.78"],
[5.54, "00:36:04.58_-10:34:33.78"],
[5.59, "00:36:02.01_-10:33:22.48"],
[6.51, "00:35:54.28_-10:32:11.18"],
[7.17, "00:35:56.86_-10:30:59.87"],
[5.56, "00:35:59.43_-10:29:48.54"],
[5.8, "00:36:02.01_-10:28:37.20"],
[5.42, "00:36:07.16_-10:28:37.20"],
[6.08, "00:36:12.31_-10:28:37.20"],
[6.34, "00:36:17.46_-10:28:37.20"],
[6.59, "00:36:20.04_-10:29:48.54"],
[5.61, "00:36:22.61_-10:30:59.87"],
[6.04, "00:36:25.19_-10:32:11.18"],
[5.61, "00:36:22.61_-10:33:22.48"],
[5.52, "00:36:20.04_-10:34:33.78"],
[5.73, "00:36:17.46_-10:35:45.06"],
[5.78, "00:36:12.31_-10:35:45.06"],
[5.33, "00:36:07.16_-10:35:45.06"],
[6.1, "00:36:02.01_-10:35:45.06"],
[5.04, "00:35:59.43_-10:34:33.78"],
[5.64, "00:35:56.86_-10:33:22.48"]]

p1_sns = []
for g in p1_grid:  
    p1_sns.append(g[0])

for pos in p1_grid:
    raj, decj = pos[1].split("_")
    coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
    dec = coord.dec.degree
    ra = coord.ra.degree
    sn = pos[0]
    
    if use_sns:
        if sn > 7.2:
            p1_ax.text(ra, dec, str(sn), color='black', fontsize=8, ha='center', va='center')
        circle = plt.Circle((ra, dec), p1r, color='r', fill=False, linewidth=1)#, alpha=sn/max(p1_sns))
    else:
        circle = plt.Circle((ra, dec), p1r, color='r', fill=False, linewidth=1)
    p1_ax.add_artist(circle)


# Phase 2 extended plot --------------------------------------------------------------------
p2e_grid_final = [[17.6, "00:36:14.70_-10:33:37.18"],
                  [8.5, "00:36:12.13_-10:33:37.18"],
                  [5.2, "00:36:13.41_-10:33:01.70"],
                  [21.5, "00:36:15.98_-10:33:01.70"],
                  [21.6, "00:36:17.26_-10:33:37.18"],
                  [13.6, "00:36:15.98_-10:34:12.65"],
                  [8.9, "00:36:13.41_-10:34:12.65"]]

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
        if sn > 7.:
            p2e_ax.text(ra, dec, str(sn), color='black', fontsize=8, ha='center', va='center')
        #p1_ax.text(ra, dec, str(sn), color='b', fontsize=8, ha='center', va='center')
        circle = plt.Circle((ra, dec), p2er, color='b', fill=False, linewidth=0.5)#, alpha=sn/max(p2e_sns))
        p2e_ax.add_artist(circle)
        #circle = plt.Circle((ra, dec), p2er, color='b', fill=False, linewidth=0.5)#, alpha=sn/max(p2e_sns))
    #else:
    #    circle = plt.Circle((ra, dec), p2er, color='b', fill=False, linewidth=0.5)
    #p1_ax.add_artist(circle)

p1_ax.set_ylim(dec_est-p2cr/2, dec_est+p2cr/2)
p1_ax.set_xlim(ra_est-p2cr/2, ra_est+p2cr/2)
#p1_ax.yaxis.tick_right()
p1_ax.set_aspect('equal', adjustable='box')

if use_sns:
    coord = SkyCoord("00:36:14.70", "-10:33:37.18", unit=(u.hourangle,u.deg))
    dec = coord.dec.degree
    ra = coord.ra.degree
    p2e_ax.set_ylim(dec-2*p2er, dec+2*p2er)
    x_ax.set_xlim(ra-2*p2er, ra+2*p2er)
    p2e_ax.set_xlim(ra-2*p2er, ra+2*p2er)
    y_ax.set_ylim(dec-2*p2er, dec+2*p2er)
    #p2e_ax.yaxis.tick_right()
    p2e_ax.set_aspect('equal', adjustable='box')

coord = SkyCoord("00:36:14.61", "-10:33:42.81", unit=(u.hourangle,u.deg))
dec = coord.dec.degree
ra = coord.ra.degree
circle = plt.Circle((ra, dec), p1r/2, color=ma.colorAlpha_to_rgb('r', 0.5)[0], zorder=0.5)
p1_ax.add_artist(circle)
circle = plt.Circle((ra, dec), p1r/2, color='r', fill=False)
p2e_ax.add_artist(circle)


hdul = fits.open('J0036-1033_post.fits')
fits_data = hdul[0].data

# Ryan size estimate
pos_u = 12 / 1.5 / 60 / 60

# Ryan position estimate
coord = SkyCoord("00:36:16.05", "-10:33:31.55", unit=(u.hourangle,u.deg))
dec = coord.dec.degree
ra = coord.ra.degree

ra_including_offset = ra - 42/60/60
dec_including_offset = dec - 11.2/60/60
#p2e_ax.arrow(ra_including_offset, dec_including_offset, 42/60/60, 11.2/60/60, head_width=0.0005, head_length=0.001, fc='r', ec='r', width=0.0001)
#print(ra_including_offset, dec_including_offset, 42/60/60, 11.2/60/60)

map_dec_range = range(30)
ra_range = np.arange(dec-2*pos_u, dec+2*pos_u-pos_u/30, 4*pos_u/30)
map_ra_range = range(30)
dec_range = np.arange(ra-2*pos_u, ra+2*pos_u-pos_u/30, 4*pos_u/30)
RA = []; Dec = []; fits_sn = []
dec_sums = np.zeros(30)
ra_sums  = np.zeros(30)
for i in map_dec_range:
    for j in map_ra_range:
        Dec.append(dec_range[i])
        RA.append(ra_range[j])
        dec_sums[i] += fits_data[i][j]
        ra_sums[j]  += fits_data[i][j]
        fits_sn.append(fits_data[i][j])
        norm_loc = round(fits_data[i][j]/np.amax(fits_data), 4)
        if norm_loc != 0:
            print("")
            print(norm_loc)
            print(ma.colorAlpha_to_rgb('b', norm_loc)[0])
        new_colour = ma.colorAlpha_to_rgb('b', norm_loc)[0]
        if new_colour[-1] > 1.:
            new_colour = np.array([new_colour[0], new_colour[1], 1.])
            print(new_colour)
        try:
            p2e_ax.scatter(dec_range[i], ra_range[j], marker='s',
                           color=ma.colorAlpha_to_rgb('b', norm_loc)[0], zorder=0.5*norm_loc)
        except:
            print("above failed")

y_ax.plot(dec_sums/max(dec_sums), np.arange(dec-2*pos_u, dec+2*pos_u-pos_u/30, 4*pos_u/30))
x_ax.plot(np.arange(ra-2*pos_u, ra+2*pos_u-pos_u/30, 4*pos_u/30), ra_sums/max(ra_sums))
x_ax.set_ylabel(r'$\rho$ (A.U.)')
y_ax.set_xlabel(r'$\rho$ (A.U.)')

p2c_ax.set_xlabel("Right Acension (degrees)")
p2c_ax.set_ylabel("Declination (degrees)")
p1_ax.set_xlabel("Right Acension (degrees)")
p2e_ax.set_xlabel("Right Acension (degrees)")

#mark zoomed in subplot
mark_inset(p2c_ax, p1_ax, loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)
mark_inset(p1_ax, p2e_ax, loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)

plt.savefig("first_mwa_discovery_localisation.png", dpi=200)
plt.savefig("first_mwa_discovery_localisation.eps")