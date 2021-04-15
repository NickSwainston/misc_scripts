import csv
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, exp, log

#with open('data/1275758864_position_localisation.csv') as csvfile:
with open('data/1275758864_position_localisation_best_pos.csv') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    data = []
    for row in spamreader:
        pulsar, status, psrcat_pointing, psrcat_prepfold_SN, psrcat_python_SN, localisation_pointing, localisation_prepfold_SN, localisation_python_SN, improvement, comment, second_comment = row
        #print(row)
        try:
            float(improvement)
            data.append(row)
        except:
            print("Skipping {} because {}".format(pulsar, comment))
        #print(', '.join(row))

positions_calcs = []
# For obsid 1275758864
fwhm = 1.05 /60. #in degrees
print("FWHM(``): {:.1f}".format(fwhm*3600.))
for cand in data:
    pulsar, status, psrcat_pointing, psrcat_prepfold_SN, psrcat_python_SN, localisation_pointing, localisation_prepfold_SN, localisation_python_SN, improvement, comment, second_comment = cand
    if float(improvement) > 0.:
        c = SkyCoord( psrcat_pointing.split("_")[0], psrcat_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
        og_ra  = c.ra.deg
        og_dec = c.dec.deg
        
        c = SkyCoord( localisation_pointing.split("_")[0], localisation_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
        loc_ra  = c.ra.deg
        loc_dec = c.dec.deg

        ra_diff = loc_ra - og_ra
        dec_diff = loc_dec - og_dec
        #print("diff(ra dec ``): {:6.2f} {:6.2f}".format(ra_diff*60.*60., dec_diff*60.*60.,))
        positions_calcs.append([pulsar, og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement])
        # work out expected improvement
        distance = np.sqrt(ra_diff**2 + dec_diff**2)
        s=fwhm/2
        x=distance*1.177
        print("")
        print("Pulsar: {}  distance(``): {:.1f}".format(pulsar, distance*3600))
        print("expected SN improvement(%): {:.2f}".format( ( (1/(exp(-(x)**2/(2*s**2)))) - 1 )*100 ))
        print("calculated  improvement(%): {}".format(improvement))
        #print(og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement)
    else:
        c = SkyCoord( psrcat_pointing.split("_")[0], psrcat_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
        og_ra  = c.ra.deg
        og_dec = c.dec.deg
        positions_calcs.append([pulsar, og_ra, og_dec, og_ra, og_dec, 0., 0., 0.])
#print(np.array(positions_calcs))

# Make a plot of the relative changes
colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
ci = 0
size = 4
plt.figure(figsize=(size, size))
for pos in positions_calcs:
    pulsar, og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement = pos
    if ra_diff > 0.:
        plt.arrow(0., 0., ra_diff*60.*60., dec_diff*60.*60., head_width=1.5, head_length=2.5, fc=colours[ci], ec=colours[ci])
        plt.scatter(60, 60, label=pulsar)
        ci += 1
pos_size = 60
plt.xlim(-pos_size, pos_size)
plt.ylim(-pos_size, pos_size)
plt.legend(loc='upper left')
plt.xlabel("Difference in Right Ascension (arcseconds)")
plt.ylabel("Difference in Declination (arcseconds)")
plt.savefig("relative_positions.png")

# Plot relative positions on the sky
size = 10
plt.figure(figsize=(size, size))
for pos in positions_calcs:
    pulsar, og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement = pos
    if ra_diff > 0.:
        scale = 100.
        #print(og_ra, og_dec, ra_diff*scale, dec_diff*scale)
        distance = np.sqrt((ra_diff*scale)**2 + (dec_diff*scale)**2)
        plt.arrow(og_ra, og_dec, ra_diff*scale, dec_diff*scale, head_width=distance*0.05, head_length=distance*0.1)#, fc=colours[ci], ec=colours[ci])
    else:
        plt.scatter(og_ra, og_dec)
plt.xlim(260., 295.)
plt.ylim( -15.,2.)
plt.xlabel("Right Ascension (degrees)")
plt.ylabel("Declination (degrees)")
#plt.legend(loc='upper left')
plt.savefig("sky_relative_positions.png", dpi=1000)