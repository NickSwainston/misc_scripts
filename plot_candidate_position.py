import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u


fig, ax = plt.subplots()

#FHRM radius
p2cr = 18.56 / 60. /2.
p1r  =  2.53 / 60. /2.

#first detection
coord = SkyCoord(['00:35:29.7', '00:36:36.23'], ['-10:35:14.03', '-10:35:14.03'],
                 unit=(u.hourangle,u.deg))
dec_centre = coord.dec.degree
ra_centre = coord.ra.degree
circle = plt.Circle((dec_centre[0], ra_centre[0]), p2cr, color='g', fill=False)
ax.add_artist(circle)

#second detection
circle = plt.Circle((dec_centre[1], ra_centre[1]), p2cr, color='g', fill=False)
ax.add_artist(circle)

#p2c posititon estimate
coord = SkyCoord("00:36:14.25", "-10:35:14.03", unit=(u.hourangle,u.deg))
dec_est = coord.dec.degree
ra_est = coord.ra.degree
circle = plt.Circle((dec_est, ra_est), 0.05, color='g', fill=False)
ax.add_artist(circle)

plt.xlim(min(dec_centre)-p2cr, max(dec_centre)+p2cr)
plt.ylim(min(ra_centre)-p2cr, max(ra_centre)+p2cr)
plt.xlabel("Right Acension (degrees)")
plt.ylabel("Declination (degrees)")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("p2c_localisation.png")
#plt.show()



plt.clf()
plt.close()
fig, ax = plt.subplots()


#zoom in
circle = plt.Circle((dec_est, ra_est), 0.05, color='g', fill=False)
ax.add_artist(circle)

p1_pos = ["00:36:14.25_-10:35:14.03","00:36:06.43_-10:35:14.03","00:36:10.34_-10:33:25.93",
          "00:36:18.15_-10:33:25.93","00:36:22.06_-10:35:14.03","00:36:18.15_-10:37:02.12",
          "00:36:10.34_-10:37:02.12"]

for pos in p1_pos:
    raj, decj = pos.split("_")
    coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
    dec = coord.dec.degree
    ra = coord.ra.degree
    circle = plt.Circle((dec, ra), p1r, color='b', fill=False)
    ax.add_artist(circle)

plt.xlim(dec_est-0.05, dec_est+0.05)
plt.ylim(ra_est-0.05, ra_est+0.05)
plt.xlabel("Right Acension (degrees)")
plt.ylabel("Declination (degrees)")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("p1_localisation.png")
#plt.show()


p1_pos2 = ["00:36:10.34_-10:33:25.93","00:36:07.89_-10:33:25.93","00:36:09.11_-10:32:52.14","00:36:11.56_-10:32:52.14","00:36:12.78_-10:33:25.93","00:36:11.56_-10:33:59.71","00:36:09.11_-10:33:59.71","00:36:05.45_-10:33:25.93","00:36:06.67_-10:32:52.14","00:36:07.89_-10:32:18.34","00:36:10.34_-10:32:18.34","00:36:12.78_-10:32:18.34","00:36:14.00_-10:32:52.14","00:36:15.22_-10:33:25.93","00:36:14.00_-10:33:59.71","00:36:12.78_-10:34:33.50","00:36:10.34_-10:34:33.50","00:36:07.89_-10:34:33.50","00:36:06.67_-10:33:59.71"]

for pos in p1_pos2:
    if pos == "00:36:14.00_-10:33:59.71":
        fill=True
    else:
        fill=False
    raj, decj = pos.split("_")
    coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
    dec = coord.dec.degree
    ra = coord.ra.degree
    circle = plt.Circle((dec, ra), p1r/4, color='b', fill=fill)
    ax.add_artist(circle)

plt.savefig("p1_localisation2.png")

coord = SkyCoord("00:36:14.61", "-10:33:42.81", unit=(u.hourangle,u.deg))
dec = coord.dec.degree
ra = coord.ra.degree
circle = plt.Circle((dec, ra), 0.0025, color='r', fill=True)
ax.add_artist(circle)
plt.savefig("p2e_localisation.png")

