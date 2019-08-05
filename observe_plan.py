from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime as DT

def time_to_observe(times, ra, dec=-26.):
    obstimes = Time(times, format='isot', scale='utc') 

    sky_posn = SkyCoord(ra, dec, unit=(u.deg,u.deg))
    #earth_location = EarthLocation.of_site('Murchison Widefield Array') 
    earth_location = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8) 
    start_Za = []
    end_Za = []
    for obstime in obstimes:
        altaz = sky_posn.transform_to(AltAz(obstime=obstime, location=earth_location)) 
        Alt = altaz.alt.deg
        if Alt[0] > 0.:
            start_Za.append(90. - Alt[0])
        else:
            start_Za.append(0.)
        if Alt[1] > 0.:
            end_Za.append(90. - Alt[1])
        else:
            end_Za.append(0.)
    print("    Date: {}".format(str(obstimes[0])[:10]))
    start_time = int(str(obstimes[start_Za.index(max(start_Za))])[11:13]) -4
    end_time   = int(str(obstimes[end_Za.index(max(end_Za))])[11:13]) -4
    
    #hours away from midnight
    start_delta_midnight = abs(start_time - 12.)
    end_delta_midnight = abs(end_time - 12.)

    am_or_pm = "pm"
    if start_time >= 12:
        am_or_pm = "am"
    if start_time > 12:
        start_time -= 12
        start_delta_midnight = start_time + \
                float(str(obstimes[start_Za.index(max(start_Za))])[14:16])/60.
    else:
        start_delta_midnight = 12. - (start_time + \
                float(str(obstimes[start_Za.index(max(start_Za))])[14:16])/60.)
        
    start_str = "{}:{} {}".format(start_time, 
                                  str(obstimes[start_Za.index(max(start_Za))])[14:16],
                                  am_or_pm)
    am_or_pm = "pm"
    if end_time >= 12:
        am_or_pm = "am"
    if end_time > 12:
        end_time -= 12
        end_delta_midnight = end_time + \
                float(str(obstimes[end_Za.index(max(end_Za))])[14:16])/60.
    else:
        end_delta_midnight = 12. - (end_time + \
                float(str(obstimes[end_Za.index(max(end_Za))])[14:16])/60.)

    end_str = "{}:{} {}".format(end_time,
                                str(obstimes[end_Za.index(max(end_Za))])[14:16],
                                am_or_pm)

    print("        start time:   {0}    end time:   {1}".format(start_str, end_str))
    print("        start zenith: {0:4.2f}      end zenith: {1:4.2f}".format(max(start_Za), max(end_Za)))
    return [start_delta_midnight, end_delta_midnight]

if __name__ == "__main__":
    #fig, ax = plt.subplots()
    #plt.gca().invert_yaxis()
    #for colour_pos in [["red",    30.,  105.], ["green",  106., 180.], ["yellow", 181., 255.]]:
    for colour_pos in [["red",   72.,  105.], ["green",  144., 180.], ["purple", 216., 255.]]:
        print("TIMES FOR {}".format(colour_pos[0]))
        year = 2019
        months = range(7, 16)
        days = range(1,30, 1)
        hours = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
        minutes = range(0, 60, 1)
        days_x_axis = np.array([])
        start_delta_midnight = []
        end_delta_midnight = []
    
        for m in months:
            if m > 12:
                year = 2020
                m -= 12
            for d in days:
                #days_x_axis.append('{0:4d}-{1:02d}-{2:02d}'.format(year, m, d))
                days_x_axis = np.append(days_x_axis, DT.datetime(year, m, d).strftime("%Y-%m-%d"))

                times = []
                for h in hours:
                    for minu in minutes:
                        times.append('{0:4d}-{1:02d}-{2:02d}T{3:02d}:{4:02d}:00'.format(year, m, d, h, minu))
                temp = time_to_observe(times, [colour_pos[1], colour_pos[2]])
                start_delta_midnight.append(temp[0])
                end_delta_midnight.append(temp[1])
        fig, ax = plt.subplots()
        plt.gca().invert_yaxis()
        days_x_axis = np.array(days_x_axis, dtype='datetime64')
        #days_x_axis = range(len(start_delta_midnight))
        """
        import matplotlib.transforms as mtransforms
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

        ax.fill_between(days_x_axis, 4, 0,
                        where = np.minimum(start_delta_midnight, end_delta_midnight) <= 2.,
                        facecolor=colour_pos[0], alpha=0.5, transform=trans)
        """
        #inbetween = np.mean([start_delta_midnight, end_delta_midnight], axis=0)
        inbetween = np.array(start_delta_midnight)
        ax.fill_between(days_x_axis, inbetween + 2.33, inbetween - 2.33, 
                         facecolor=colour_pos[0], alpha=0.5)
        ax.plot(days_x_axis, inbetween + 2.33, colour_pos[0])
        ax.plot(days_x_axis, inbetween - 2.33,
                colour_pos[0], label="RA: {} deg".format(colour_pos[1]))
        """
        ax.plot(days_x_axis, start_delta_midnight, colour_pos[0],
                label="RA: {} deg".format(colour_pos[1]))
        if colour_pos[0] is "red":
            colour_alt = "darkred"
        elif colour_pos[0] is "green":
            colour_alt = "darkgreen"
        elif colour_pos[0] is "purple":
            colour_alt = "indigo"
        ax.plot(days_x_axis, end_delta_midnight, colour_alt,
                label="RA: {} deg".format(colour_pos[2]))
        """
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        #ax.xaxis.set_major_locator(mdates.DayLocator())
        plt.setp(ax.get_xticklabels(), rotation=45)
        plt.gcf().autofmt_xdate()
        ax.set_ylim([4.,0.])
        plt.legend()
        plt.xlabel('Date')
        plt.ylabel('Hours away from midnight')
        #plt.show()
        plt.savefig("{}_group.png".format(colour_pos[0]))
        plt.clf()
