import numpy as np

import matplotlib.dates as mdates
from astroplan import AtNightConstraint, AltitudeConstraint
from astroplan.utils import time_grid_from_range
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astroplan.plots import plot_airmass
from astroplan import Observer, FixedTarget
from astropy.coordinates import SkyCoord
from matplotlib import colors
import datetime as DT

ra = 180.
dec = -26.
red_ra_decs    = [[43.4,1.6],[43.4,-26.7],[50.7,-55.0],[59.8,-40.5],[59.8,18.3],[59.8,-13.0],[70.7,-72.0], [76.2,1.6],[76.2,-26.7],[90.8,-55.0],[92.6,-40.5],[92.6,18.3],[92.6,-13.0]]
green_ra_decs  = [[109.0,1.6],[109.0,-26.7],[125.4,-40.5],[125.4,18.3],[125.4,-13.0],[130.9,-55.0],[130.9,-72.0], [141.8,1.6],[141.8,-26.7],[158.2,-40.5],[158.2,18.3],[158.2,-13.0],[171.0,-55.0],[174.6,1.6], [174.6,-26.7]]
#purple_ra_decs = [[190.1,-13.0],[190.1,18.3],[190.1,-72.0],[190.1,-40.5],[206.5,-26.7],[206.5,1.6],[210.1,-55.0], [222.9,-13.0],[222.9,18.3],[222.9,-40.5],[239.3,-26.7],[239.3,1.6],[250.2,-72.0],[250.2,-55.0]]
purple_ra_decs = [[239.3,-26.7],[239.3,1.6],[250.2,-72.0],[250.2,-55.0]]
yellow_ra_decs = [[255.7,-13.0],[255.7,18.3],[255.7,-40.5],[272.1,-26.7],[272.1,1.6],[288.5,-13.0],[288.5,18.3], [288.5,-40.5],[290.3,-55.0],[304.9,-26.7],[304.9,1.6],[310.4,-72.0],[321.3,-13.0],[321.3,18.3],[321.3,-40.5]]
blue_ra_decs   = [[330.4,-55.0], [337.7,1.6],[337.7,-26.7],[354.1,-13.0],[354.1,18.3],[354.1,-40.5],[10.5,-26.7], [10.5,1.6],[10.6,-55.0],[10.6,-72.0],[26.9,-40.5],[26.9,18.3],[26.9,-13.0]]

# Set up observation targets
colour_targets = []
for ci, colour_ra_dec in enumerate([["purple", purple_ra_decs], ["orange", yellow_ra_decs]]): #["red", red_ra_decs], ["green", green_ra_decs], ["blue", blue_ra_decs]
    #targets = []
    colour, ra_dec = colour_ra_dec
    for ra, dec in ra_dec:
        target = FixedTarget(SkyCoord( ra, dec, unit=(u.deg,u.deg)), name=colour)
        colour_targets.append([colour, target])
    #colour_targets.append(targets)
observer = Observer(longitude=116.670813889*u.deg, latitude=-26.703319444*u.deg, elevation=377.8*u.m)

# Set up time range to check
start_time = Time('2021-11-01T16:00:01', format='isot', scale='utc')
end_time   = Time('2022-12-31T16:00:01', format='isot', scale='utc')
time_resolution = 24 * u.hour
time_grid = time_grid_from_range([start_time, end_time],
                                 time_resolution=time_resolution)


# Set pair of constraints
alt_cut_off1 = 20.
constraint1 = AltitudeConstraint(min=alt_cut_off1 * u.deg)
alt_cut_off2 = 40.
constraint2 = AltitudeConstraint(min=alt_cut_off2* u.deg)

# Loop over targets, check them and plot the results
plot_data = np.zeros((20, len(time_grid))) + 1
days_x_axis = np.array([])
fig, ax = plt.subplots(figsize=(15,5))
for ci, colour_target in enumerate(colour_targets):
    colour, target = colour_target

    # loose constraint
    obs1 = constraint1(observer, target, times=time_grid)
    valid_times = time_grid[obs1]
    plt.scatter(min(valid_times.value), ci, alpha=0)
    rectangle = plt.Rectangle((min(valid_times.value), ci+0.5), max(valid_times.value) - min(valid_times.value), 1, fc='white', ec=colour)
    plt.gca().add_patch(rectangle)

    # tight constraint
    obs2 = constraint2(observer, target, times=time_grid)
    valid_times = time_grid[obs2]
    plt.scatter(min(valid_times.value), ci, alpha=0)
    rectangle = plt.Rectangle((min(valid_times.value), ci+0.5), max(valid_times.value) - min(valid_times.value), 1, fc=colour)
    plt.gca().add_patch(rectangle)

ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.xlabel("Date (month-day)")
#plt.setp(ax.get_xticklabels(), rotation=45)
#plt.gcf().autofmt_xdate()
#fig.autofmt_xdate()
#plt.yticks([], [])
ax.set_yticks(list(range(1, len(colour_targets))))
plt.title("Date ranges the SMART observations are above an altitude of {} deg for empty rectangles and {} deg for filled rectangles".format(alt_cut_off1, alt_cut_off2))
plt.savefig('SMART_observing_plan.png')


