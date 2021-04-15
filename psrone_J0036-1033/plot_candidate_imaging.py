import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
from astropy.visualization import astropy_mpl_style
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
#plt.style.use(astropy_mpl_style)

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord


cmap='gray_r'
labelsize = 6
#plt.rcParams["font.family"] = "Times New Roman"
#fontname = 'Times New Roman'
plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False
fig, (p2c_ax, p1_ax, p2e_ax) = plt.subplots(1,3, figsize=(9,3))
#fig = plt.figure(figsize=(11, 3))

#index, position, centre, size (arcsec), [min Jy, max Jy], file
data_array = [#[1, ['00h36m15.15', '−10d33m11.90'], ['0h36m16.05', '-10d33m31.55'], 3600, [-0.08, 0.22],  '/astro/mwavcs/rbhat/psrone/fits/downsampled_avg_1278106408_1275085816_I_NEW2.fits'],
              #[2, ['00h36m15.15', '−10d33m11.90'], ['0h36m16.05', '-10d33m31.55'], 600,  [-1e-5, 15e-5], '/astro/mwavcs/rbhat/psrone/fits/J0036-1022_band3.fits'],
              #[3, ['00h36m14.97', '-10d33m16.02'], ['0h36m16.05', '-10d33m31.55'], 300,  [-5e-6, 24e-6], '/astro/mwavcs/rbhat/psrone/fits/J0036-1033_03dec20_band4.fits']]
              #[3, ['00h36m14.97', '-10d33m16.02'], ['0h36m16.05', '-10d33m31.55'], 300,  [-5e-6, 24e-6], '/astro/mwavcs/rbhat/psrone/fits/J0036-1033_650MHz_1.fits']]
              [1, ['00h36m15.15', '−10d33m11.90'], ['0h36m16.05', '-10d33m31.55'], 3600, [-0.08, 0.22],  '/astro/mwavcs/nswainston/code/misc_scripts/psrone_J0036-1033/J0036-1033_MWA_si.ci.image_new.fits'],
              [2, ['00h36m15.15', '−10d33m11.90'], ['0h36m16.05', '-10d33m31.55'], 600,  [-1e-4, 20e-4], '/astro/mwavcs/nswainston/code/misc_scripts/psrone_J0036-1033/J0036-1033-GMRTB3.ci.image.fits'],
              [3, ['00h36m14.97', '-10d33m16.02'], ['0h36m16.05', '-10d33m31.55'], 275,  [-5e-5, 24e-5], '/astro/mwavcs/nswainston/code/misc_scripts/psrone_J0036-1033/J0036-1033-GMRTB4.ci.image.fits']]
axes = []
wcss = []
box_pixs = []
for index, position, centre, size, jy, file_loc in data_array:
    centre = position
    hdu = fits.open(file_loc)[0]
    print("index {}".format(index))
    print("hdu.data.shape {}".format(len(hdu.data.shape)))
    #for h in hdu.header:
    #    print(h, hdu.header[h])
    wcs = WCS(hdu.header, naxis=2)
    size = u.Quantity((size, size), u.arcsec)
    if len(hdu.data.shape) == 4:
        s_data = np.squeeze(hdu.data)
    else:
        s_data = hdu.data
    print("s_data.shape {}".format(s_data.shape))
    data = Cutout2D(s_data, SkyCoord(position[0], position[1], frame='icrs'), size, fill_value=np.nan, wcs=wcs)
    #print(data)
    ax = plt.subplot(1, 3, index, projection=data.wcs)
    #ax = WCSAxes(fig, [0.1+index-1, 0.1, 0.8, 0.8], wcs=data.wcs)
    #fig.add_axes(ax)  
    axes.append(ax)
    wcss.append(data.wcs)
    #if index == 1:
    #    data_factor = 1000
    #else:
    #    data_factor = 10000
    data_factor = 1000
    im = ax.imshow(data.data*data_factor, cmap=cmap,
                   aspect='equal',
                   origin='lower', 
                   interpolation='none',
                   rasterized=True,
                   vmax=jy[1]*data_factor,
                   vmin=jy[0]*data_factor)
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=None)
    cbar.ax.minorticks_on()
    cbar.update_ticks()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(labelsize)
    cbar.set_label('mJy/beam', size=labelsize)
    ax.set_xlabel("Right Acension (J2000)", size=labelsize)
    ax.set_ylabel("Declination (J2000)", size=labelsize)
    lon, lat = ax.coords
    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)

    # add cand pocition circle
    circle_size = data.shape[0] * 0.04
    pixel_position = data.wcs.world_to_pixel(SkyCoord(position[0], position[1], frame='icrs'))
    circle = plt.Circle(pixel_position, circle_size, color='r', lw=0.5, fill=False)
    ax.add_patch(circle)

    # make a box
    if index > 1:
        print("index {} skycord: {}".format(index, SkyCoord(centre[0], centre[1], frame='icrs')))
        centre_position = data.wcs.world_to_pixel(SkyCoord(centre[0], centre[1], frame='icrs'))
        print("index {} centre_position: {}".format(index, centre_position))
        left  = centre_position[0] - data.shape[0]/2
        right = centre_position[0] + data.shape[0]/2
        top   = centre_position[1] + data.shape[1]/2
        bot   = centre_position[1] - data.shape[1]/2
        # TL TR BR BL TL
        box_x = [left, right, right, left, left]
        box_y = [top, top, bot, bot, top]
        box_pixs.append([box_x, box_y])

    dec = ax.coords[1]
    dec.set_ticklabel(size=6)
    ra = ax.coords[0]
    ra.set_ticklabel(size=6)

#for ax in axes:
#        labels = ax.get_xticklabels() + ax.get_yticklabels()
#        [label.set_fontname(fontname) for label in labels]
#mark zoomed in subplot
#position = SkyCoord(data_array[0][1][0], data_array[0][1][1], frame='icrs')
# Create a set of inset Axes: these should fill the bounding box allocated to them.
#axes[1] = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
#ip = InsetPosition(axes[0], [0.4,0.2,0.5,0.5])
#axes[1].set_axes_locator(ip)
#mark_inset(axes[0], axes[1], loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)#, transform=wcs[0])
#indicate_inset(axes[0], [], fc="none", ec="0.5", zorder=0.5)#, loc1=2, loc2=3, transform=wcs[0])
#mark_inset(axes[1], axes[2], loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)

# manually mark box
box_lines = []
for i, box in enumerate(box_pixs):
    # convert to ra dec from pixels of previos wcs
    box_pix_x = []
    box_pix_y = []
    for pos in zip(box[0], box[1]):
        #print(wcss[i+1])
        sky = wcss[i+1].pixel_to_world(pos[0], pos[1])
        x, y = wcss[i].world_to_pixel(sky)
        box_pix_x.append(x)
        box_pix_y.append(y)
    box_lines.append([box_pix_x, box_pix_y])
print(box_lines)
axes[0].plot(box_lines[0][0], box_lines[0][1], lw=0.5)
axes[1].plot(box_lines[1][0], box_lines[1][1], lw=0.5)

#fig.tight_layout()
plt.subplots_adjust(wspace=0.8)
plt.savefig("plot_candidate_imaging.png", dpi=1200)
plt.savefig("plot_candidate_imaging.eps", format='eps', dpi=500)
