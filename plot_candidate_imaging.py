import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
from astropy.visualization import astropy_mpl_style
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.style.use(astropy_mpl_style)

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord


cmap='gray_r'
labelsize = 8
fig, (p2c_ax, p1_ax, p2e_ax) = plt.subplots(1,3, figsize=(11,3))
#fig = plt.figure(figsize=(11, 3))

#index, position, size, file
data_array = [[1, ['0h36m20.0', '-10d35m00.00'], 3600, '/astro/mwavcs/rbhat/psrone/fits/downsampled_avg_1278106408_1275085816_I_NEW2.fits'],
              [2, ['0h36m15.0', '-10d33m00.00'], 600, '/astro/mwavcs/rbhat/psrone/fits/J0036-1022_band3.fits'],
              [3, ['0h36m15.0', '-10d33m30.00'], 300, '/astro/mwavcs/rbhat/psrone/fits/J0036-1033_03dec20_band4.fits']]

axes = []
wcss = []
for index, position, size, file_loc in data_array:
    hdu = fits.open(file_loc)[0]
    #for h in hdu.header:
    #    print(h, hdu.header[h])
    wcs = WCS(hdu.header, naxis=2)
    size = u.Quantity((size, size), u.arcsec) 
    data = Cutout2D(hdu.data, SkyCoord(position[0], position[1], frame='icrs'), size, fill_value=np.nan, wcs=wcs)
    print(data.shape)
    ax = plt.subplot(1, 3, index, projection=data.wcs)
    #ax = WCSAxes(fig, [0.1+index-1, 0.1, 0.8, 0.8], wcs=data.wcs)
    #fig.add_axes(ax)  
    axes.append(ax)
    wcss.append(data.wcs)
    im = ax.imshow(data.data, cmap=cmap,
                   aspect = 'equal',
                   origin = 'lower', 
                   interpolation = 'none',
                   rasterized = True)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(8)
    cbar.set_label('Jy/beam', size=labelsize)
    ax.set_xlabel("Right Acension (J2000)", size=labelsize)
    ax.set_ylabel("Declination (J2000)", size=labelsize)

    dec = ax.coords[1]
    dec.set_ticklabel(size=6)
    ra = ax.coords[0]
    ra.set_ticklabel(size=6)

#mark zoomed in subplot
position = SkyCoord(data_array[0][1][0], data_array[0][1][1], frame='icrs')
# Create a set of inset Axes: these should fill the bounding box allocated to them.
#axes[1] = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
#ip = InsetPosition(axes[0], [0.4,0.2,0.5,0.5])
#axes[1].set_axes_locator(ip)
mark_inset(axes[0], axes[1], loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)#, transform=wcs[0])
#indicate_inset(axes[0], [], fc="none", ec="0.5", zorder=0.5)#, loc1=2, loc2=3, transform=wcs[0])
mark_inset(axes[1], axes[2], loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)

#fig.tight_layout()
plt.subplots_adjust(wspace=0.8)
plt.savefig("plot_candidate_imaging.png", dpi=200)