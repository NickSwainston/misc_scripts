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
red_target    = FixedTarget(SkyCoord( 67., dec, unit=(u.deg,u.deg)), name="red")
green_target  = FixedTarget(SkyCoord(139., dec, unit=(u.deg,u.deg)), name="green")
purple_target = FixedTarget(SkyCoord(211., dec, unit=(u.deg,u.deg)), name="purple")
yellow_target = FixedTarget(SkyCoord(283., dec, unit=(u.deg,u.deg)), name="yellow")
blue_target   = FixedTarget(SkyCoord(354., dec, unit=(u.deg,u.deg)), name="blue")
observer = Observer(longitude=116.670813889*u.deg, latitude=-26.703319444*u.deg, elevation=377.8*u.m)

start_time = Time('2019-07-01T16:00:01', format='isot', scale='utc')
end_time   = Time('2020-08-01T16:00:01', format='isot', scale='utc')
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


