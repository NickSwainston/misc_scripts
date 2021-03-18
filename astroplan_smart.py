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
purple_ra_decs = [[190.1,-13.0],[190.1,18.3],[190.1,-72.0],[190.1,-40.5],[206.5,-26.7],[206.5,1.6],[210.1,-55.0], [222.9,-13.0],[222.9,18.3],[222.9,-40.5],[239.3,-26.7],[239.3,1.6],[250.2,-72.0],[250.2,-55.0]]
yellow_ra_decs = [[255.7,-13.0],[255.7,18.3],[255.7,-40.5],[272.1,-26.7],[272.1,1.6],[288.5,-13.0],[288.5,18.3], [288.5,-40.5],[290.3,-55.0],[304.9,-26.7],[304.9,1.6],[310.4,-72.0],[321.3,-13.0],[321.3,18.3],[321.3,-40.5]]
blue_ra_decs   = [[330.4,-55.0], [337.7,1.6],[337.7,-26.7],[354.1,-13.0],[354.1,18.3],[354.1,-40.5],[10.5,-26.7], [10.5,1.6],[10.6,-55.0],[10.6,-72.0],[26.9,-40.5],[26.9,18.3],[26.9,-13.0]]

colour_targets = []
for ci, colour_ra_dec in enumerate([["red", red_ra_decs], ["green", green_ra_decs], ["purple", purple_ra_decs], ["yellow", yellow__ra_decs], ["blue", blue_ra_decs]]):
    targets = []
    colour, ra_dec = colour_ra_dec
    for ra, dec in ra_dec:
        target = FixedTarget(SkyCoord( ra, dec, unit=(u.deg,u.deg)), name="red")
        targets.append(target)
    colour_targets.append(targets)
#red_target    = FixedTarget(SkyCoord( 67., dec, unit=(u.deg,u.deg)), name="red")
#green_target  = FixedTarget(SkyCoord(139., dec, unit=(u.deg,u.deg)), name="green")
#purple_target = FixedTarget(SkyCoord(211., dec, unit=(u.deg,u.deg)), name="purple")
#yellow_target = FixedTarget(SkyCoord(283., dec, unit=(u.deg,u.deg)), name="yellow")
#blue_target   = FixedTarget(SkyCoord(354., dec, unit=(u.deg,u.deg)), name="blue")
observer = Observer(longitude=116.670813889*u.deg, latitude=-26.703319444*u.deg, elevation=377.8*u.m)

start_time = Time('2021-05-01T16:00:01', format='isot', scale='utc')
end_time   = Time('2021-12-31T16:00:01', format='isot', scale='utc')
time_resolution = 24 * u.hour

time_grid = time_grid_from_range([start_time, end_time],
                                 time_resolution=time_resolution)

#red_nights = AtNightConstraint(observer, target, times=time_grid,
#                                         max_solar_altitude=30. * u.deg)
#red_alts = AltitudeConstraint(observer, target, times=time_grid,
#                                        min=60. * u.deg)
#red_obs = np.logical_and(red_nights, red_alts)
constraints = [AltitudeConstraint(min=50. * u.deg)]
for i, constraint in enumerate(constraints):
    red_obs    = constraint(observer, red_target,    times=time_grid)
    green_obs  = constraint(observer, green_target,  times=time_grid)
    purple_obs = constraint(observer, purple_target, times=time_grid)
    yellow_obs = constraint(observer, yellow_target, times=time_grid)
    blue_obs   = constraint(observer, blue_target,   times=time_grid)

plot_data = np.zeros((5, len(time_grid))) + 1
days_x_axis = np.array([])
for t, time in enumerate(time_grid):
    if red_obs[t]:
        plot_data[0, t] = 3
    if green_obs[t]:
        plot_data[1, t] = 5
    if purple_obs[t]:
        plot_data[2, t] = 7
    if yellow_obs[t]:
        plot_data[3, t] = 9
    if blue_obs[t]:
        plot_data[4, t] = 11
    days_x_axis = np.append(days_x_axis, time.datetime.strftime("%Y-%m-%d %H:%M:%S"))
days_x_axis = np.array(days_x_axis, dtype='datetime64')

cmap = colors.ListedColormap(['white','red', 'green', 'purple', 'yellow', 'blue'])
bounds = [0,2,4,6,8, 10, 12]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(15,5))
ax.imshow(plot_data, extent=[mdates.date2num(days_x_axis[0]), mdates.date2num(days_x_axis[-1]), 0, 2],
                             aspect='auto', cmap=cmap, norm=norm)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
#plt.setp(ax.get_xticklabels(), rotation=45)
#plt.gcf().autofmt_xdate()
fig.autofmt_xdate()
plt.yticks([], [])
# draw gridlines
#ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
xticks = []
for d in days_x_axis:
    if "01T16:00:00" in str(d):
        xticks.append(d)
ax.set_xticks(xticks);
#ax.set_yticks(np.arange(-.5, 10, 1));

plt.savefig('SMART_observing_plan.png')


