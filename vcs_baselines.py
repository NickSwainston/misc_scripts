
from astropy.io import fits
import numpy as np
import seaborn
import matplotlib.pyplot as plt

def getTileLocations(metafits, flags=[]):
    """
    Function grab the MWA tile locations for any given observation ID. Downloads the relevant metafits file from the database, saves it as <obdID>_metafits_ppds.fits.

    Input:
      obsid - the GPS observation ID
      flags - RTS tile flags (i.e. the first entry in the metafits correspond to "tile 0", irrespective of what the antenna name is)

    Return:
      a list of lists containing the following:
        list[0] = a list of tile positions East of the array centre
        list[1] = a list of tile positions North of the array centre
        list[2] = a list of tile heights about sea-level
    """

    f = fits.open(metafits)

    east = f[1].data['East'][::2]
    north = f[1].data['North'][::2]
    height = f[1].data['Height'][::2] # height above sea-level

    # MWA array centre height above sea-level
    mwacentre_h = 377.827
    height = height - mwacentre_h

    # flag the tiles from the x,y,z positions
    east = np.delete(east,flags)
    north = np.delete(north,flags)
    height = np.delete(height,flags)

    return east, north, height

if __name__ == "__main__":
    # phase 2 compact example
    metafits = '/astro/mwavcs/vcs/1223042480/1223042480_metafits_ppds.fits'
    east, north, height = getTileLocations(metafits, flags=[])
    print(east)

    # loop over the tiles to get all the baselines
    baselines = []
    for i in range(len(east)):
        for j in range(i+1, len(east)):
            # Calculate distance
            dist = np.sqrt( (east[i] - east[j])**2 + (north[i] - north[j])**2 + (height[i] - height[j])**2 )
            baselines.append(dist)

    #print(baselines)
    print(np.mean(baselines))
    print(np.median(baselines))

    # count baselines less than 360
    count = 0
    for b in baselines:
        if b < 360.:
            count += 1
    print(count)
    print(count/len(baselines)*100)

    baseline_plot = seaborn.displot(baselines)
    #fig = baseline_plot.get_figure()
    plt.savefig("baseline_distribution.png") 
