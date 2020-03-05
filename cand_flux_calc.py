import os
import argparse
import numpy as np
import subprocess
import sys
from shutil import copyfile as cp
import math
import glob
import textwrap as _textwrap
from astropy.table import Table

#MWA software imports
from mwa_pulsar_client import client
from mwa_metadb_utils import get_common_obs_metadata
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe
import prof_utils
import sys

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

trcvr="/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"
bestprof_data = prof_utils.get_from_bestprof(sys.argv[1])
obsid, pulsar, _, period, _, beg, t_int, profile, num_bins = bestprof_data
period=float(period)
num_bins=int(num_bins)

#get r_sys and gain
obs_metadata = get_common_obs_metadata(obsid)
p_ra = "00:36:14.25"
p_dec = "10:35:14.03"

obsid, obs_ra, obs_dec, _, delays, centrefreq, channels = obs_metadata

#get beg if not supplied
if beg is None or t_int is None:
    logger.debug("Calculating beginning time for pulsar coverage")
    beg, _, t_int = find_times(obsid, pulsar, beg=beg, end=end, min_z_power=min_z_power)

#Find 'start_time' for fpio - it's usually about 7 seconds
#obs_start, _ = mwa_metadb_utils.obs_max_min(obsid)
start_time = beg-int(obsid)

#Get important info
trec_table = Table.read(trcvr,format="csv")
ntiles = 128 #TODO actually we excluded some tiles during beamforming, so we'll need to account for that here

beam_power = fpio.get_beam_power_over_time([obsid, obs_ra, obs_dec, t_int, delays,\
                                            centrefreq, channels],\
                                            np.array([[pulsar, p_ra, p_dec]]),\
                                            dt=100, start_time=start_time)
beam_power = np.mean(beam_power)

# Usa a primary beam function to convolve the sky temperature with the primary beam
# (prints suppressed)
sys.stdout = open(os.devnull, 'w')
from mwa_pb import primarybeammap_tant as pbtant
_, _, Tsky_XX, _, _, _, Tsky_YY, _ = pbtant.make_primarybeammap(int(obsid), delays, centrefreq*1e6, 'analytic', plottype='None')
sys.stdout = sys.__stdout__


#TODO can be inaccurate for coherent but is too difficult to simulate
t_sky = (Tsky_XX + Tsky_YY) / 2.
import submit_to_database
import mwa_metadb_utils
# Get T_sys by adding Trec and Tsky (other temperatures are assumed to be negligible
t_sys_table = t_sky + submit_to_database.get_Trec(trec_table, centrefreq)
t_sys = np.mean(t_sys_table)
t_sys_err = t_sys*0.02 #TODO: figure out what t_sys error is

logger.debug("pul_ra: {} pul_dec: {}".format(p_ra, p_dec))
_, _, zas = mwa_metadb_utils.mwa_alt_az_za(obsid, ra=p_ra, dec=p_dec)
theta = np.radians(zas)
gain = submit_to_database.from_power_to_gain(beam_power, centrefreq*1e6, ntiles, coh=True)
logger.debug("beam_power: {} theta: {} pi: {}".format(beam_power, theta, np.pi))
gain_err = gain * ((1. - beam_power)*0.12 + 2.*(theta/(0.5*np.pi))**2. + 0.1)

# Removed the below error catch because couldn't find an obs that breaks it
#sometimes gain_err is a numpy array and sometimes it isnt so i have to to this...
#try:
#    gain_err.shape
#    gain_err = gain_err[0]
#estimate S/N
try:
    prof_dict = prof_utils.auto_gfit(profile,\
                period = period, plot_name="{0}_{1}_{2}_bins_gaussian_fit.png".format(obsid, pulsar, num_bins))
except prof_utils.ProfileLengthError:
    prof_dict=None

if not prof_dict:
    logger.info("Profile couldn't be fit. Using old style of profile analysis")
    sn, u_sn, _, w_equiv_bins, u_w_equiv_bins, w_equiv_ms, u_w_equiv_ms, scattering, u_scattering, scattered =\
        prof_utils.analyse_pulse_prof(profile, period, verbose=True)
    w_equiv_ms = period/num_bins * w_equiv_bins
    u_w_equiv_ms = period/num_bins * u_w_equiv_bins
else:
    sn = prof_dict["sn"]
    u_sn = prof_dict["sn_e"]
    w_equiv_bins = prof_dict["Weq"]
    u_w_equiv_bins =  prof_dict["Weq_e"]
    w_equiv_ms = period/num_bins * w_equiv_bins
    u_w_equiv_ms = period/num_bins * u_w_equiv_bins
    scattering = prof_dict["Wscat"]*period/num_bins/1000 #convert to seconds
    u_scattering = prof_dict["Wscat_e"]*period/num_bins/1000
    scattered = prof_dict["scattered"]

logger.info("Profile scattered? {0}".format(scattered))
logger.info("S/N: {0} +/- {1}".format(sn, u_sn))
logger.debug("Gain {0} K/Jy".format(gain))
logger.debug("Equivalent width in bins: {0}".format(w_equiv_bins))
logger.debug("T_sys: {0} K".format(t_sys))
bandwidth = 30720000. #In Hz
logger.debug("Bandwidth: {0}".format(30.72))
time_detection = t_int
logger.debug("Detection time: {0}".format(time_detection))
logger.debug("NUmber of bins: {0}".format(num_bins))

if scattered == False:
    #final calc of the mean fluxdesnity in mJy
    S_mean = sn * t_sys / ( gain * math.sqrt(2. * float(time_detection) * bandwidth)) *\
             math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
    #constants to make uncertainty calc easier
    S_mean_cons = t_sys / ( math.sqrt(2. * float(time_detection) * bandwidth)) *\
             math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
    u_S_mean = math.sqrt( math.pow(S_mean_cons * u_sn / gain , 2)  +\
                          math.pow(sn * S_mean_cons * gain*0.1 / math.pow(gain,2) , 2) )

    logger.info('Smean {0:.2f} +/- {1:.2f} mJy'.format(S_mean, u_S_mean))
