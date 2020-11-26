from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, ICRS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import math
import numpy as np

from mwa_metadb_utils import get_common_obs_metadata, get_obs_array_phase, getmeta
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe
import prof_utils

import logging
logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)
#ch = logging.StreamHandler()
#ch.setLevel(logging.WARNING)

detections = [ [1255444104, "/astro/mwavcs/nswainston/J0036-1033_detections/1255444104_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 41.55],
          [1133775752, "/astro/mwavcs/nswainston/J0036-1033_detections/1133775752_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 8.38],
          [1137236608, "/astro/mwavcs/nswainston/J0036-1033_detections/1194350120_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 10.54],
          [1150234552, "/astro/mwavcs/nswainston/J0036-1033_detections/1150234552_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 5.44],
          [1164110416, "/astro/mwavcs/nswainston/J0036-1033_detections/1164110416_00:36:14.25_-10:35:14.03_900.00ms_Cand.pfd.bestprof", 36.31],
          [1182630616, "/astro/mwavcs/nswainston/J0036-1033_detections/1182630616_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 6.4],
          [1194350120, "/astro/mwavcs/nswainston/J0036-1033_detections/1194350120_00:36:14.61_-10:33:42.81_900.00ms_Cand.pfd.bestprof", 8.01],
          [1220886016, "/astro/mwavcs/nswainston/J0036-1033_detections/1220886016_00:36:14.00_-10:33:59.71_ch109-120_900.00ms_Cand.pfd.bestprof", 15.36],
          [1220886016, "/astro/mwavcs/nswainston/J0036-1033_detections/1220886016_00:36:14.00_-10:33:59.71_ch145-156_900.00ms_Cand.pfd.bestprof", 12.91],
          [1222697776, "/astro/mwavcs/nswainston/J0036-1033_detections/1222697776_00:36:10.34_-10:33:25.93_900.00ms_Cand.pfd.bestprof", 31.53],
          [1225462936, "/astro/mwavcs/nswainston/J0036-1033_detections/1225462936_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 12.8],
          [1275085816, "/astro/mwavcs/nswainston/J0036-1033_detections/1275085816_00:36:11.58_-10:33:56.44_900.04ms_Cand.pfd.bestprof", 17.92],
          [1275092416, "/astro/mwavcs/nswainston/J0036-1033_detections/1275092416_00:36:12.34_-10:32:48.36_900.04ms_Cand.pfd.bestprof", 17.85],
          [1275094456, "/astro/mwavcs/nswainston/J0036-1033_detections/1275094456_00:36:13.85_-10:33:26.18_900.04ms_Cand.pfd.bestprof", 14.62],
          [1275172216, "/astro/mwavcs/nswainston/J0036-1033_detections/1275172216_00:36:15.87_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 20.78],
          [1275177136, "/astro/mwavcs/nswainston/J0036-1033_detections/1275177136_00:36:14.78_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 19.65],
          [1275178816, "/astro/mwavcs/nswainston/J0036-1033_detections/1275178816_00:36:15.20_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 16.94],
          [1275258616, "/astro/mwavcs/nswainston/J0036-1033_detections/1275258616_00:36:14.61_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 17.66],
          [1275431416, "/astro/mwavcs/nswainston/J0036-1033_detections/1275431416_00:36:09.67_-10:33:15.91_900.04ms_Cand.pfd.bestprof", 24.08],
          [1275863416, "/astro/mwavcs/nswainston/J0036-1033_detections/1275863416_00:36:12.76_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 25.52],
          [1275866536, "/astro/mwavcs/nswainston/J0036-1033_detections/1275866536_00:36:13.68_-10:33:24.92_900.04ms_Cand.pfd.bestprof", 28.5],
          [1276725752, "/astro/mwavcs/nswainston/J0036-1033_detections/1276725752_00:36:13.60_-10:33:48.87_900.04ms_Cand.pfd.bestprof", 18.53],
          [1278106408, "/astro/mwavcs/nswainston/J0036-1033_detections/1278106408_00:36:14.44_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 27.36],
          [1283104232, "/astro/mwavcs/nswainston/J0036-1033_detections/1283104232_00:36:14.02_-10:33:22.44_900.04ms_Cand.pfd.bestprof", 24.07],
          [1285086000, "/astro/mwavcs/nswainston/J0036-1033_detections/1285086000_00:36:13.77_-10:33:22.44_900.04ms_Cand.pfd.bestprof", 17.84],
          [1287670032, "/astro/mwavcs/nswainston/J0036-1033_detections/1287670032_00:36:14.02_-10:33:33.78_900.04ms_Cand.pfd.bestprof", 25.52]]


pulsar="J0036-1033"
bestprof_data = prof_utils.get_from_bestprof("/astro/mwavcs/nswainston/J0036-1033_detections/1278106408_00:36:14.44_-10:33:18.61_900.04ms_Cand.pfd.bestprof")
obsid, prof_psr, dm, period, _, beg, t_int, profile, num_bins = bestprof_data
prof_dict = prof_utils.auto_gfit(profile, period)


mjds = [58774.602847222224, 57366.404340277775, 57406.460543981484, 57556.89971064815, 57717.49998842592, 57931.85414351852, 58067.49655092593, 58374.624976851854, 58374.624976851854, 58395.59442129629, 58427.598587962966, 59001.937476851854, 59002.013865740744, 59002.03747685185, 59002.937476851854, 59002.994421296295, 59003.013865740744, 59003.937476851854, 59005.937476851854, 59010.937476851854, 59010.973587962966, 59020.918217592596, 59036.89803240741, 59094.74321759259, 59117.680347222224, 59147.588125]
sns = [36.180054684322116, 8.611114663132401, 11.451918495336324, 10.725666576249761, 20.165451103011378, 8.673979646933207, 6.668467172336903, 11.256089731789622, 11.18038554673332, 21.638526831026724, 13.152361686289439, 15.527820290146945, 15.26124674940845, 10.423596979620521, 13.672227372811143, 14.95603184468743, 12.099388406216525, 12.155173188105056, 17.15820357106891, 17.137561143963175, 19.325804317013546, 12.741559707702663, 20.5133490037488, 15.49284802740998, 13.27150743247289, 16.930259903441325]
norm_sns = [11.185885277272176, 6.608036692245963, 3.938236272714118, 1.9688263602442928, 3.347324769739635, 4.748671794035306, 3.3165869299190986, 5.2913753944942705, 7.123091432736747, 8.057288501091563, 4.088049056095448, 5.854663701572106, 8.047772834209637, 4.733012412119289, 5.337611741157288, 5.787081219398446, 5.417902834811325, 6.178442309783125, 5.679017331670005, 6.980214700524004, 7.641481362864146, 4.201340912437257, 5.859816077210512, 4.944149595113389, 4.52659570319818, 4.3010402585646315]
u_norm_sns = [0.09301068692306166, 0.8738270534163127, 0.3332504550652208, 0.20368716153703356, 0.10965035289864414, 0.5660108017894304, 0.6909981869378609, 0.4693127751707024, 0.6682025683674854, 0.20429417247275294, 0.2943755488126769, 0.27493411489815833, 0.35524617099203665, 0.38163252135348846, 0.31219613024867277, 0.25338221422951773, 0.31681043574435425, 0.35116499504868925, 0.18587513643431863, 0.2321007109369484, 0.19609987672565432, 0.24123151332659423, 0.1481330050375039, 0.19614663157240814, 0.22898043791152783, 0.14318285725370777]
colours = ['purple', 'r', 'r', 'r', 'g', 'g', 'b', 'g', 'g', 'g', 'g', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b']
array_phases = ['P2C', 'P1', 'P1', 'P1', 'P2C', 'P2C', 'P2E', 'P2C', 'P2C', 'P2C', 'P2C', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E']

'''
array_phases = []
colours = []
mjds = []
sns = []
norm_sns = []
u_norm_sns = []
pdmp_norm_sns = []
print("obsid, sn, u_sn_factional, sn_normalised, u_sn_normalised")
for det in detections:
    obsid, bestprof_file, pdmp_sn = det
    mjd = Time(int(obsid), format='gps', scale='utc').mjd
    o_phase = get_obs_array_phase(obsid)
    if o_phase == 'OTH':
        o_phase = 'P2E'
    if obsid == 1255444104:
        #orig obs
        colours.append('purple')
    elif o_phase == 'P2E':
        colours.append('b')
    elif o_phase == 'P2C':
        colours.append('g')
    elif o_phase == 'P1':
        colours.append('r')
    array_phases.append(o_phase)
    #print(obsid)
    #unpack the bestprof_data
    #[obsid, pulsar, dm, period, period_uncer, obsstart, obslength, profile, bin_num]
    bestprof_data = prof_utils.get_from_bestprof(bestprof_file)
    obsid, prof_psr, dm, period, _, beg, t_int, profile, num_bins = bestprof_data
    if not pulsar:
        pulsar = prof_psr
    period=float(period)
    num_bins=int(num_bins)

    #get r_sys and gain
    metadata = get_common_obs_metadata(obsid)
    _, ra_obs, dec_obs, time_obs, delays, centrefreq, channels = metadata
    if bestprof_file == "/astro/mwavcs/nswainston/J0036-1033_detections/1275092416_00:36:12.34_-10:32:48.36_900.04ms_Cand.pfd.bestprof":
        minfreq = 139.52
        maxfreq = 153.60
        bandwidth = 15360000. #In Hz
    elif  bestprof_file ==  "/astro/mwavcs/nswainston/J0036-1033_detections/1220886016_00:36:14.00_-10:33:59.71_ch145-156_900.00ms_Cand.pfd.bestprof":
        minfreq = 185.60
        maxfreq = 199.68
        bandwidth = 15360000. #In Hz
    else:
        minfreq = float(min(channels)) * 1.28
        maxfreq = float(max(channels)) * 1.28
        bandwidth = 30720000. #In Hz

    #Get minimum beam offset
    beam_distance = []
    cand_posn = SkyCoord("00:36:14.70", "-10:33:37.18", unit=(u.hourangle, u.deg))
    earth_location = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
    beam_meta_data = getmeta(service='obs', params={'obs_id':obsid})
    alt = beam_meta_data[u'metadata'][u'elevation_pointing'] #in sexidecimal
    az = beam_meta_data[u'metadata'][u'azimuth_pointing']
    #beam_atlaz = SkyCoord(az, alt, unit=(u.hourangle, u.deg), frame='altaz')#, unit=(u.alt.deg,u.az.deg))
    #beam_atlaz = SkyCoord(alt, az, unit=(u.hourangle, u.deg), frame='altaz')
    for obstime_raw in range(obsid, obsid + time_obs, 100):
        obstime = Time(float(obstime_raw),format='gps')
        #altaz = cand_posn.transform_to(AltAz(obstime=obstime, location=earth_location))
        beam_atlaz = SkyCoord(az, alt, unit=(u.deg, u.deg), frame='altaz',
                              obstime=obstime, location=earth_location)
        tile_radec = beam_atlaz.transform_to(ICRS())
        #print(tile_radec.ra.deg, tile_radec.dec.deg)
        #print(tile_radec.separation(cand_posn).deg)
        #beam_distance.append(beam_atlaz.separation(altaz).deg)
        beam_distance.append(tile_radec.separation(cand_posn).deg)
    min_beam_offset = min(beam_distance)
    
    output_data, obsid_meta = fpio.find_sources_in_obs([obsid], [["J0036-1033", "00:36:14.70", "-10:33:37.18"]],
                                 dt_input=100, min_power=0.3, degrees_check=False)
    pulsar, enter, exit, max_power = output_data[obsid][0]
    #pdmp_sn_normalised = pdmp_sn / ( max_power * math.sqrt(float(t_int)/4800 * bandwidth/30720000)) 
    #normalise for 20 mins
    pdmp_sn_normalised = pdmp_sn / ( max_power * math.sqrt(float(t_int)/1200 * bandwidth/30720000))


    #t_sys, _, gain, u_gain = snfe.find_t_sys_gain(pulsar, obsid, obs_metadata=metadata,\
    #                                p_ra="00:36:14.70", p_dec="-10:33:37.18",
    #                                beg=beg, end=(t_int + beg - 1))

    #estimate S/N
    """
    try:
        prof_dict = prof_utils.auto_gfit(profile,\
                    period = period, plot_name="{0}_{1}_{2}_bins_gaussian_fit.png".format(obsid, pulsar, num_bins))
    except (prof_utils.ProfileLengthError, prof_utils.NoFitError) as _:
        prof_dict=None
    sn = prof_dict["sn"]
    u_sn = prof_dict["sn_e"]
    """
    #Normalize, find the SN
    # Normalise
    y = np.array(profile)/max(profile)
    # Find std of the noise
    noise_std, clipped = prof_utils.sigmaClip(y, alpha=3.)

    # re normalise it fo 0 is on the noise mean
    y = y - np.nanmean(clipped)
    y = y/max(y)    
    # Simple method
    """
    # Sum on pulse
    pulse_flux = 0.
    pulse_width = 0
    for yi in range(len(y)):
        #print(clipped[yi])
        if np.isnan(clipped[yi]):
            pulse_flux = pulse_flux + y[yi]
            pulse_width = pulse_width + 1
            #print(pulse_flux)
    print(pulse_width)
    """
    # Make sure we get all the on pulse method
    max_index = list(y).index(max(y))
    pulse_width = 10
    rotate_by = len(y) // 2 - max_index
    y = list(y[rotate_by:]) + list(y[:rotate_by])
    max_index = list(y).index(max(y))
    pulse_flux = sum(y[max_index-5:max_index+5])

    # centre pulse
    if pulse_flux == 0:
        print(y[max_index-5:max_index+5])
        print(y, max_index)
    sn = pulse_flux / pulse_width / noise_std
    sn_peak = max(y) / noise_std
    u_sn_factional = np.sqrt(pulse_width) * noise_std / sn 
    #sn = pulse_flux / period / noise_std
    #u_sn = np.sqrt(period) * noise_std / sn

    sn_normalised = sn / ( max_power * math.sqrt(float(t_int)/1200 * bandwidth/30720000))
    u_sn_normalised = sn_normalised * u_sn_factional

    print("{} {:5.1f} {:5.1f} {:5.1f} {:5.1f}".format(obsid, sn, u_sn_factional, sn_normalised, u_sn_normalised))

    """
    #sn = pdmp_sn
    #u_sn = pdmp_sn * 0.1
    w_equiv_bins = prof_dict["Weq"]
    u_w_equiv_bins =  prof_dict["Weq_e"]
    w_equiv_ms = period/num_bins * w_equiv_bins
    u_w_equiv_ms = period/num_bins * u_w_equiv_bins
    scattering = prof_dict["Wscat"]*period/num_bins/1000 #convert to seconds
    u_scattering = prof_dict["Wscat_e"]*period/num_bins/1000
    scattered = prof_dict["scattered"]

    
    #if scattered:
    #    print(f'{obsid} & {mjd} & scattered')
    #    continue
    #print("S/N: {0} +/- {1}".format(sn, u_sn))
    logger.debug("Gain {0} K/Jy".format(gain))
    logger.debug("Equivalent width in bins: {0}".format(w_equiv_bins))
    logger.debug("T_sys: {0} K".format(t_sys))
    logger.debug("Bandwidth: {0}".format(bandwidth))
    logger.debug("Detection time: {0}".format(t_int))
    logger.debug("NUmber of bins: {0}".format(num_bins))

    #final calc of the mean fluxdesnity in mJy
    S_mean = sn * t_sys / ( gain * math.sqrt(2. * float(t_int) * bandwidth)) *\
            math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
    #constants to make uncertainty calc easier
    S_mean_cons = t_sys / ( math.sqrt(2. * float(t_int) * bandwidth)) *\
            math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
    u_S_mean = math.sqrt( math.pow(S_mean_cons * u_sn / gain , 2)  +\
                        math.pow(sn * S_mean_cons * u_gain / math.pow(gain,2) , 2) )
    """
    
    #if obsid == 1220886016:
    #    S_mean = S_mean/2

    mjds.append(mjd)
    sns.append(sn_peak)
    norm_sns.append(sn_normalised)
    u_norm_sns.append(u_sn_normalised)
    pdmp_norm_sns.append(pdmp_sn_normalised)
    t_int = int(t_int)
    #print(f'{obsid} & {mjd:.1f} & {o_phase:3}     & {minfreq:.2f}-{maxfreq:.2f} & {t_int:4}      & {min_beam_offset:6.1f}     & {sn:5.1f} & {sn_normalised:6.1f} \\\\')
    #{S_mean:.2f} & {u_S_mean:.2f} \\\\') 
    #print('Smean {0:.2f} +/- {1:.2f} mJy'.format(S_mean, u_S_mean))
'''


print("mjds = {}".format(mjds))
print("sns = {}".format(sns))
print("norm_sns = {}".format(norm_sns))
print("u_norm_sns = {}".format(u_norm_sns))
print("colours = {}".format(colours))
print("array_phases = {}".format(array_phases))

print("obsid, MJD, SN, norm_SN")
for i in range(len(mjds)):
    print("{}, {:.1f}, {:5.1f}, {:5.1f}".format(detections[i][0], mjds[i], sns[i], norm_sns[i]))



# Cut off change scale -----------------
fig ,(ax1,ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', figsize=(10, 5))

#plot
#ax1.scatter(mjds[:11], norm_sns[:11])
#ax2.scatter(mjds[11:], norm_sns[11:])

markersize = 3
makerwidth = 1
capsize = 3

# plot ax1
array_phase_legend = {"P1": True, "P2C": True, "P2E": True}
array_phase_legend_labels = {"P1": "Phase 1 array", "P2C": "Phase 2 Compact Array", "P2E":"Phase 2 Extended Array"}
for i in range(11):
    if norm_sns[i] == max(norm_sns):
        (plotline, caps, barlinecols) = ax1.errorbar(mjds[i], norm_sns[i], yerr=u_norm_sns[i], c=colours[i], fmt="o", markersize=markersize, capsize=capsize, label="1st Pulsar Detection")
    else:
        #print(array_phase_legend[array_phases[i]])
        (plotline, caps, barlinecols) = ax1.errorbar(mjds[i], norm_sns[i], yerr=u_norm_sns[i], c=colours[i], fmt="o", markersize=markersize, capsize=capsize,
                                                     label=array_phase_legend_labels[array_phases[i]] if array_phase_legend[array_phases[i]] else "")
        array_phase_legend[array_phases[i]] = False
    for cap in caps:
        cap.set_markeredgewidth(makerwidth)
    if norm_sns[i] == max(norm_sns):
        #orig obs
        od = (plotline, caps)
    elif array_phases[i] == "P1":
        p1pl = (plotline, caps)
    elif array_phases[i] == "P2C":
        p2cpl = (plotline, caps)
    elif array_phases[i] == "P2E":
        p2epl = (plotline, caps)
for i in range(11, len(mjds)):
    #if i in [22, 11]:
    #    colours[i] = 'orange'
    #    print(colours[i])
    (_, caps, _) = ax2.errorbar(mjds[i], norm_sns[i], yerr=u_norm_sns[i], c=colours[i], fmt="o", markersize=markersize, capsize=capsize)
    for cap in caps:
        cap.set_markeredgewidth(makerwidth)
    #ax2.legend(array_phases[i])
#ax2.legend((od, p1pl, p2cpl, p2epl), ('1st', 'P1', 'P2C', 'P2E'))
ax1.legend(loc='upper left')
# hide the spines between ax and ax2
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
#ax1.tick_params(labelright='off')
ax2.yaxis.set_label_position("right")
ax2.tick_params(axis='y', which='both', labelleft='off', labelright='on')
ax2.yaxis.tick_right()

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d,1+d), (-d,+d), **kwargs)
ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d), (1-d,1+d), **kwargs)
ax2.plot((-d,+d), (-d,+d), **kwargs)

# Set up tics
"""
ax1.set_xticks(np.arange(57000, 59000, 10), minor = True)
ax1.set_xticks(np.arange(57000, 59000, 100))
ax2.set_xticks(np.arange(59000, 59040, 10), minor = True)
ax2.set_xticks(np.arange(59000, 59040, 100))
print(np.arange(57000, 59000, 10))
print(np.arange(57000, 59000, 100))
print(np.arange(59000, 59040, 10))
print(np.arange(59000, 59040, 100))
"""
ax1.xaxis.set_major_locator(MultipleLocator(200))
ax1.xaxis.set_minor_locator(MultipleLocator(20))
ax2.xaxis.set_major_locator(MultipleLocator(50))
ax2.xaxis.set_minor_locator(MultipleLocator(5))

ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=4)


ax1.set_xticklabels((ax1.get_xticks()-57300).astype(int))
ax2.set_xticklabels((ax2.get_xticks()-59000).astype(int))
ax1.set_ylabel('Mean Flux (Arbitrary Units)')
ax1.set_xlabel('Days since MJD {}'.format(57300))
ax2.set_xlabel('Days since MJD {}'.format(59000))

plt.savefig('normalised_sn_scale_change.png', bbox_inches='tight', dpi=1000)