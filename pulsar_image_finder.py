import sys
import numpy as np

from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u

from vcstools.catalogue_utils import get_psrcat_ra_dec
from vcstools.analyse_psf import read_psf_fits


fits_file = sys.argv[1]
print("Load fits file {}".format(fits_file))
rav, decv, fits_data = read_psf_fits(fits_file)
pixel_size = decv[1][0] - decv[0][0]
rav = rav.flatten()
decv = decv.flatten()
fits_data = fits_data.flatten()
pixel_coords = SkyCoord(ra=rav, dec=decv, unit=(u.hourangle,u.deg))

print("Get pulsar ra and decs")
jname_ra_dec = np.array(get_psrcat_ra_dec())
jname = jname_ra_dec[:,0]
ra    = jname_ra_dec[:,1]
dec   = jname_ra_dec[:,2]
pulsar_coords = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))

print("Finding the pixels near a pulsar")
# Check if they're within a beam width
pulsar_idx, pixel_idx, sep2d, _ = search_around_sky(pulsar_coords, pixel_coords, pixel_size*2*u.deg)

# Loop over pulsars to check each one
for pi, pulsar in enumerate(jname):
    if pi in pulsar_idx:
        # Pulsar is in field of view so collect all pixels
        seps = []
        pidxs = []
        for pu_idx, pi_idx, sep in zip(pulsar_idx, pixel_idx, sep2d):
            if pu_idx == pi:
                seps.append(sep)
                pidxs.append(pi_idx)
        best_index = pidxs[seps.index(min(seps))]
        best_pixel = pixel_coords[best_index]
        print("{:11s} value: {:.5f}   Pixel location: {:6.2f} {:6.2f}".format(pulsar, fits_data[best_index], best_pixel.ra.deg, best_pixel.dec.deg))