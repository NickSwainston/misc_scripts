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

detections = [[1133775752, "/group/mwavcs/vcs/1133775752/pointings/00:36:14.00_-10:33:59.71/1133775752_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 8.38],
          [1137236608, "/group/mwavcs/vcs/1137236608/pointings/00:36:14.00_-10:33:59.71/1194350120_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 10.54],
          [1150234552, "/group/mwavcs/vcs/1150234552/pointings/00:36:14.00_-10:33:59.71/1150234552_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 5.44],
          [1164110416, "/group/mwavcs/vcs/1164110416/pointings/00:36:14.25_-10:35:14.03/1164110416_00:36:14.25_-10:35:14.03_900.00ms_Cand.pfd.bestprof", 36.31],
          [1182630616, "/group/mwavcs/vcs/1182630616/pointings/00:36:14.00_-10:33:59.71/1182630616_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 6.4],
          [1194350120, "/group/mwavcs/vcs/1194350120/pointings/00:36:14.61_-10:33:42.81/1194350120_00:36:14.61_-10:33:42.81_900.00ms_Cand.pfd.bestprof", 8.01],
          [1220886016, "/group/mwavcs/vcs/1220886016/pointings/00:36:14.00_-10:33:59.71/1220886016_00:36:14.00_-10:33:59.71_ch109-120_900.00ms_Cand.pfd.bestprof", 15.36],
          [1220886016, "/group/mwavcs/vcs/1220886016/pointings/00:36:14.00_-10:33:59.71/1220886016_00:36:14.00_-10:33:59.71_ch145-156_900.00ms_Cand.pfd.bestprof", 12.91],
          [1222697776, "/group/mwavcs/vcs/1222697776/pointings/00:36:10.34_-10:33:25.93/1222697776_00:36:10.34_-10:33:25.93_900.00ms_Cand.pfd.bestprof", 31.53],
          [1225462936, "/group/mwavcs/vcs/1225462936/pointings/00:36:14.00_-10:33:59.71/1225462936_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 12.8],
          [1255444104, "/group/mwavcs/vcs/1255444104/pointings/00:36:14.00_-10:33:59.71/1255444104_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 41.55],
          [1275085816, "/group/mwavcs/nswainston/pulsar_search/1275085816_candidate_follow_up/1275085816_00:36:11.58_-10:33:56.44_900.04ms_Cand.pfd.bestprof", 17.92],
          [1275092416, "/group/mwavcs/nswainston/pulsar_search/1275092416_candidate_follow_up/1275092416_00:36:12.34_-10:32:48.36_900.04ms_Cand.pfd.bestprof", 17.85],
          [1275094456, "/group/mwavcs/nswainston/pulsar_search/1275094456_candidate_follow_up/1275094456_00:36:13.85_-10:33:26.18_900.04ms_Cand.pfd.bestprof", 14.62],
          [1275172216, "/group/mwavcs/nswainston/pulsar_search/1275172216_candidate_follow_up/1275172216_00:36:15.87_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 20.78],
          [1275177136, "/group/mwavcs/nswainston/pulsar_search/1275177136_candidate_follow_up/1275177136_00:36:14.78_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 19.65],
          [1275178816, "/group/mwavcs/nswainston/pulsar_search/1275178816_candidate_follow_up/1275178816_00:36:15.20_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 16.94],
          [1275258616, "/group/mwavcs/nswainston/pulsar_search/1275258616_candidate_follow_up/1275258616_00:36:14.61_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 17.66],
          [1275431416, "/group/mwavcs/nswainston/pulsar_search/1275431416_candidate_follow_up/1275431416_00:36:09.67_-10:33:15.91_900.04ms_Cand.pfd.bestprof", 24.08],
          [1275863416, "/group/mwavcs/nswainston/pulsar_search/1275863416_candidate_follow_up/1275863416_00:36:12.76_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 25.52],
          [1275866536, "/group/mwavcs/nswainston/pulsar_search/1275866536_candidate_follow_up/1275866536_00:36:13.68_-10:33:24.92_900.04ms_Cand.pfd.bestprof", 28.5],
          [1276725752, "/group/mwavcs/nswainston/pulsar_search/1276725752_candidate_follow_up/1276725752_00:36:13.60_-10:33:48.87_900.04ms_Cand.pfd.bestprof", 18.53],
          [1278106408, "/group/mwavcs/nswainston/pulsar_search/1278106408_candidate_follow_up/1278106408_00:36:14.44_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 27.36]]


pulsar="J0036-1033"
bestprof_data = prof_utils.get_from_bestprof("/group/mwavcs/nswainston/pulsar_search/1278106408_candidate_follow_up/1278106408_00:36:14.44_-10:33:18.61_900.04ms_Cand.pfd.bestprof")
obsid, prof_psr, dm, period, _, beg, t_int, profile, num_bins = bestprof_data
prof_dict = prof_utils.auto_gfit(profile, period)

mjds = []
norm_sns = []
u_norm_sns = []
pdmp_norm_sns = []
for det in detections:
    obsid, bestprof_file, pdmp_sn = det
    mjd = Time(int(obsid), format='gps', scale='utc').mjd
    o_phase = get_obs_array_phase(obsid)
    if o_phase == 'OTH':
        o_phase = 'P2E'
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
    if bestprof_file == "/group/mwavcs/nswainston/pulsar_search/1275092416_candidate_follow_up/1275092416_00:36:12.34_-10:32:48.36_900.04ms_Cand.pfd.bestprof":
        minfreq = 139.52
        maxfreq = 153.60
        bandwidth = 15360000. #In Hz
    elif  bestprof_file ==  "/group/mwavcs/vcs/1220886016/pointings/00:36:14.00_-10:33:59.71/1220886016_00:36:14.00_-10:33:59.71_ch145-156_900.00ms_Cand.pfd.bestprof":
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
    # Sum on pulse
    pulse_flux = 0.
    pulse_width = 0
    for yi in range(len(y)):
        #print(clipped[yi])
        if np.isnan(clipped[yi]):
            pulse_flux = pulse_flux + y[yi]
            pulse_width = pulse_width + 1
            #print(pulse_flux)
    sn = pulse_flux / pulse_width / noise_std
    u_sn = np.sqrt(pulse_width) * noise_std

    sn_normalised = sn / ( max_power * math.sqrt(float(t_int)/1200 * bandwidth/30720000))
    u_sn_normalised = u_sn / ( max_power * math.sqrt(float(t_int)/1200 * bandwidth/30720000))

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
    norm_sns.append(sn_normalised)
    u_norm_sns.append(u_sn_normalised)
    pdmp_norm_sns.append(pdmp_sn_normalised)
    t_int = int(t_int)
    print(f'{obsid} & {mjd:.1f} & {o_phase:3}     & {minfreq:.2f}-{maxfreq:.2f} & {t_int:4}      & {min_beam_offset:6.1f}     & {sn:5.1f} & {sn_normalised:6.1f} \\\\')
    #{S_mean:.2f} & {u_S_mean:.2f} \\\\') 
    #print('Smean {0:.2f} +/- {1:.2f} mJy'.format(S_mean, u_S_mean))

#mjds = [57366.404340277775, 57406.460543981484, 57556.89971064815, 57717.49998842592, 57931.85414351852, 58067.49655092593, 58374.624976851854, 58374.624976851854, 58395.59442129629, 58427.598587962966, 58774.602847222224, 59001.937476851854, 59002.013865740744, 59002.03747685185, 59002.937476851854, 59002.994421296295, 59003.013865740744, 59003.937476851854, 59005.937476851854, 59010.937476851854, 59010.973587962966, 59020.918217592596, 59036.89803240741]
#norm_sns = [39.880341880158696, 25.440146393870187, 7.516520454649676, 50.77783967013483, 19.872526637138165, 23.34347123156014, 51.31816682731386, 60.99878134934648, 88.15310930886803, 31.343352198895765, 88.43070646255248, 48.384597128356425, 61.20510339706202, 36.782424844791855, 56.097347505717366, 47.0961123200065, 41.06825962295899, 47.67525931070088, 48.57128515657029, 64.20016640202387, 68.31086018782896, 36.02152126125464, 52.5816839829118]


print(mjds)
print(norm_sns)
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].set_ylim(0, 100)
ax[1].set_ylim(0, 100)
ax[0].errorbar(mjds, norm_sns, yerr=u_norm_sns, fmt=".")
ax[1].errorbar(mjds[11:], norm_sns[11:], yerr=u_norm_sns[11:], fmt=".")
ax[1].set_yticks([])

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax[0], ax[1], loc1=2, loc2=3, fc="none", ec="0.5", zorder=0.5)

ax[0].set_xlabel('MJD')
ax[1].set_xlabel('MJD')
ax[0].set_ylabel('Normalised S/N')
#plt.show()
plt.savefig('normalised_sn_zoom.png', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(5, 5))
ax.set_ylim(0, 100)
ax.errorbar(mjds, norm_sns, yerr=u_norm_sns, fmt=".")
ax.set_xlabel('MJD')
ax.set_ylabel('Normalised S/N')
plt.savefig('normalised_sn.png', bbox_inches='tight')

#inset plots
fig, ax = plt.subplots()#figsize=(5, 5))

ax.set_ylim(0, 100)
ax.errorbar(mjds, norm_sns, yerr=u_norm_sns, fmt=".")
ax.set_xlabel('MJD')
ax.set_ylabel('Normalised S/N')

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Insert 1 -----------------------------------
zoom_factor = 1
axins1 = inset_axes(ax, 2, 1, loc='lower left', bbox_to_anchor=(1, 1), bbox_transform=ax.figure.transFigure)
axins1.errorbar(mjds[11:], norm_sns[11:], yerr=u_norm_sns[11:], fmt=".")
axins1.set_xlim(min(mjds[11:])-100, max(mjds[11:])+100) # apply the x-limits
axins1.set_ylim(min(norm_sns[11:])-2, max(norm_sns[11:])+2) # apply the y-limits
#axins1.set_yticks([y1[0]-0.00002, y1[0]-0.00001, y1[0], y1[0]+0.00001, y1[0]+0.00002])
#axins1.set_yticklabels(['-0.00002', '-0.00001', '0', '0.00001', '0.00002'])

axins1.set_xticks([], minor=True)

plt.savefig('normalised_sn_inset.png', bbox_inches='tight')

# Cut off change scale -----------------
fig ,(ax1,ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', figsize=(10, 5))

#plot
#ax1.scatter(mjds[:11], norm_sns[:11])
#ax2.scatter(mjds[11:], norm_sns[11:])
ax1.errorbar(mjds[:11], norm_sns[:11], yerr=u_norm_sns[:11], fmt=".")
ax2.errorbar(mjds[11:], norm_sns[11:], yerr=u_norm_sns[11:], fmt=".")

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
ax2.xaxis.set_major_locator(MultipleLocator(200))
ax2.xaxis.set_minor_locator(MultipleLocator(20))

ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=4)


plt.xlabel('MJD')
ax1.set_ylabel('Normalised S/N')

plt.savefig('normalised_sn_scale_change.png', bbox_inches='tight')