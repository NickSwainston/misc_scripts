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
import csv
import pandas as pd

from vcstools.metadb_utils import get_common_obs_metadata, get_obs_array_phase, getmeta
from vcstools.beam_calc import find_sources_in_obs
# import vcstools.sn_flux_est as snfe
import vcstools.prof_utils as prof_utils

import logging
logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)
#ch = logging.StreamHandler()
#ch.setLevel(logging.WARNING)

detections = [
    (1255444104, "/astro/mwavcs/nswainston/J0036-1033_detections/1255444104_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 41.55),
    (1133775752, "/astro/mwavcs/nswainston/J0036-1033_detections/1133775752_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 8.38),
    (1137236608, "/astro/mwavcs/nswainston/J0036-1033_detections/1194350120_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 10.54),
    (1150234552, "/astro/mwavcs/nswainston/J0036-1033_detections/1150234552_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 5.44),
    (1164110416, "/astro/mwavcs/nswainston/J0036-1033_detections/1164110416_00:36:14.25_-10:35:14.03_900.00ms_Cand.pfd.bestprof", 36.31),
    (1182630616, "/astro/mwavcs/nswainston/J0036-1033_detections/1182630616_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 6.4),
    (1194350120, "/astro/mwavcs/nswainston/J0036-1033_detections/1194350120_00:36:14.61_-10:33:42.81_900.00ms_Cand.pfd.bestprof", 8.01),
    # (1220886016, "/astro/mwavcs/nswainston/J0036-1033_detections/1220886016_00:36:14.00_-10:33:59.71_ch109-120_900.00ms_Cand.pfd.bestprof", 15.36),
    # (1220886016, "/astro/mwavcs/nswainston/J0036-1033_detections/1220886016_00:36:14.00_-10:33:59.71_ch145-156_900.00ms_Cand.pfd.bestprof", 12.91),
    (1222697776, "/astro/mwavcs/nswainston/J0036-1033_detections/1222697776_00:36:10.34_-10:33:25.93_900.00ms_Cand.pfd.bestprof", 31.53),
    (1225462936, "/astro/mwavcs/nswainston/J0036-1033_detections/1225462936_00:36:14.00_-10:33:59.71_900.00ms_Cand.pfd.bestprof", 12.8),
    (1275085816, "/astro/mwavcs/nswainston/J0036-1033_detections/1275085816_00:36:11.58_-10:33:56.44_900.04ms_Cand.pfd.bestprof", 17.92),
    (1275092416, "/astro/mwavcs/nswainston/J0036-1033_detections/1275092416_00:36:12.34_-10:32:48.36_900.04ms_Cand.pfd.bestprof", 17.85),
    (1275094456, "/astro/mwavcs/nswainston/J0036-1033_detections/1275094456_00:36:13.85_-10:33:26.18_900.04ms_Cand.pfd.bestprof", 14.62),
    (1275172216, "/astro/mwavcs/nswainston/J0036-1033_detections/1275172216_00:36:15.87_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 20.78),
    (1275177136, "/astro/mwavcs/nswainston/J0036-1033_detections/1275177136_00:36:14.78_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 19.65),
    (1275178816, "/astro/mwavcs/nswainston/J0036-1033_detections/1275178816_00:36:15.20_-10:33:31.22_900.04ms_Cand.pfd.bestprof", 16.94),
    (1275258616, "/astro/mwavcs/nswainston/J0036-1033_detections/1275258616_00:36:14.61_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 17.66),
    (1275431416, "/astro/mwavcs/nswainston/J0036-1033_detections/1275431416_00:36:09.67_-10:33:15.91_900.04ms_Cand.pfd.bestprof", 24.08),
    (1275863416, "/astro/mwavcs/nswainston/J0036-1033_detections/1275863416_00:36:12.76_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 25.52),
    (1275866536, "/astro/mwavcs/nswainston/J0036-1033_detections/1275866536_00:36:13.68_-10:33:24.92_900.04ms_Cand.pfd.bestprof", 28.5),
    (1276725752, "/astro/mwavcs/nswainston/J0036-1033_detections/1276725752_00:36:13.60_-10:33:48.87_900.04ms_Cand.pfd.bestprof", 18.53),
    (1278106408, "/astro/mwavcs/nswainston/J0036-1033_detections/1278106408_00:36:14.44_-10:33:18.61_900.04ms_Cand.pfd.bestprof", 27.36),
    (1283104232, "/astro/mwavcs/nswainston/J0036-1033_detections/1283104232_00:36:14.02_-10:33:22.44_900.04ms_Cand.pfd.bestprof", 24.07),
    (1285086000, "/astro/mwavcs/nswainston/J0036-1033_detections/1285086000_00:36:13.77_-10:33:22.44_900.04ms_Cand.pfd.bestprof", 17.84),
    (1287670032, "/astro/mwavcs/nswainston/J0036-1033_detections/1287670032_00:36:14.02_-10:33:33.78_900.04ms_Cand.pfd.bestprof", 25.52),
    (1290341112, "/astro/mwavcs/nswainston/J0036-1033_detections/1290341112_00:36:14.58_-10:33:16.40_900.04ms_Cand.pfd.bestprof", 30.23),
]

df = pd.DataFrame(
    columns=[
        "obsid",
        "mjd",
        "o_phase",
        "freq (MHz)",
        "flux",
        "flux_error",
        "bestprof",
        "pdmp_sn",
    ]
)

for obsid, bestprof_file, pdmp_sn in detections:
    with open(f"/astro/mwavcs/nswainston/J0036-1033_detections/None_{obsid}_dummy_flux_results.csv") as f:
        reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONE)
        for row in reader:
            # print(row)
            if row[0] == "flux":
                flux = float(row[1])
            if row[0] == "flux_error":
                flux_error = float(row[1])
    mjd = Time(int(obsid), format='gps', scale='utc').mjd
    freq = get_common_obs_metadata(obsid)[5]
    o_phase = get_obs_array_phase(obsid)
    if o_phase == 'OTH':
        o_phase = 'P2E'
    if obsid == 1255444104:
        o_phase = 'orig'

    df = df.append(
        {
            "obsid": obsid,
            "mjd": mjd,
            "o_phase": o_phase,
            "freq (MHz)": freq,
            "flux": flux,
            "flux_error": flux_error,
            "bestprof": bestprof_file,
            "pdmp_sn": pdmp_sn,
        },
        ignore_index = True
    )

print(df)
df.to_csv('psrone_variability.csv')

# mjds = [58774.602847222224, 57366.404340277775, 57406.460543981484, 57556.89971064815, 57717.49998842592, 57931.85414351852, 58067.49655092593, 58374.624976851854, 58374.624976851854, 58395.59442129629, 58427.598587962966, 59001.937476851854, 59002.013865740744, 59002.03747685185, 59002.937476851854, 59002.994421296295, 59003.013865740744, 59003.937476851854, 59005.937476851854, 59010.937476851854, 59010.973587962966, 59020.918217592596, 59036.89803240741, 59094.74321759259, 59117.680347222224, 59147.588125]
# sns = [36.180054684322116, 8.611114663132401, 11.451918495336324, 10.725666576249761, 20.165451103011378, 8.673979646933207, 6.668467172336903, 11.256089731789622, 11.18038554673332, 21.638526831026724, 13.152361686289439, 15.527820290146945, 15.26124674940845, 10.423596979620521, 13.672227372811143, 14.95603184468743, 12.099388406216525, 12.155173188105056, 17.15820357106891, 17.137561143963175, 19.325804317013546, 12.741559707702663, 20.5133490037488, 15.49284802740998, 13.27150743247289, 16.930259903441325]
# norm_sns = [11.185885277272176, 6.608036692245963, 3.938236272714118, 1.9688263602442928, 3.347324769739635, 4.748671794035306, 3.3165869299190986, 5.2913753944942705, 7.123091432736747, 8.057288501091563, 4.088049056095448, 5.854663701572106, 8.047772834209637, 4.733012412119289, 5.337611741157288, 5.787081219398446, 5.417902834811325, 6.178442309783125, 5.679017331670005, 6.980214700524004, 7.641481362864146, 4.201340912437257, 5.859816077210512, 4.944149595113389, 4.52659570319818, 4.3010402585646315]
# u_norm_sns = [0.09301068692306166, 0.8738270534163127, 0.3332504550652208, 0.20368716153703356, 0.10965035289864414, 0.5660108017894304, 0.6909981869378609, 0.4693127751707024, 0.6682025683674854, 0.20429417247275294, 0.2943755488126769, 0.27493411489815833, 0.35524617099203665, 0.38163252135348846, 0.31219613024867277, 0.25338221422951773, 0.31681043574435425, 0.35116499504868925, 0.18587513643431863, 0.2321007109369484, 0.19609987672565432, 0.24123151332659423, 0.1481330050375039, 0.19614663157240814, 0.22898043791152783, 0.14318285725370777]
# colours = ['purple', 'r', 'r', 'r', 'g', 'g', 'b', 'g', 'g', 'g', 'g', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b']
# marks = ['x', 'o', 'o', 'o', 's', 's', '^', 's', 's', 's', 's', '^', '^', '^', '^', '^', '^', '^', '^', '^', '^', '^', '^', '^', '^', '^', '^']
# array_phases = ['P2C', 'P1', 'P1', 'P1', 'P2C', 'P2C', 'P2E', 'P2C', 'P2C', 'P2C', 'P2C', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E', 'P2E']



# Cut off change scale -----------------
fig ,(ax1, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', figsize=(10, 5))


markersize = 5
marker_border_thickness = 0.5
makerwidth = 1
capsize = 3
errorbar_linewidth = 0.7

# plot ax1
array_phase_legend = {"P1": True, "P2C": True, "P2E": True}
array_phase_legend_labels = {"P1": "Phase 1 Array", "P2C": "Phase 2 Compact Array", "P2E":"Phase 2 Extended Array"}
for index, df_row in df.iterrows():
    if df_row["obsid"] < 1260000000:
        # archive
        ax = ax1
    else:
        ax = ax2

    # Set up label for first instance
    o_phase = df_row["o_phase"]
    if o_phase == 'orig':
        label = "1st Pulsar Detection"
    else:
        if array_phase_legend[o_phase]:
            label = array_phase_legend_labels[o_phase]
        else:
            label = ""
        array_phase_legend[o_phase] = False

    # Set up markers
    if o_phase == 'orig':
        #orig obs
        colour = 'purple'
        mark = "X"
    elif o_phase == 'P2E':
        colour = 'b'
        mark = "^"
    elif o_phase == 'P2C':
        colour = 'g'
        mark = "H"
    elif o_phase == 'P1':
        colour = 'r'
        mark = "o"

    # plot
    print(
        df_row["mjd"],
        df_row["flux"],
        df_row["flux_error"],
        colour,
        mark,
        markersize,
        capsize,
        label,
    )
    (plotline, caps, barlinecols) = ax.errorbar(
        df_row["mjd"],
        df_row["flux"],
        yerr=df_row["flux_error"],
        c=colour,
        fmt=mark,
        markersize=markersize,
        markeredgewidth=marker_border_thickness,
        capsize=capsize,
        elinewidth=errorbar_linewidth,
        label=label,
    )
    for cap in caps:
        cap.set_markeredgewidth(makerwidth)

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
ax1.set_ylabel('Flux Density (mJy)')
ax1.set_xlabel('Days since MJD {} (2015 Oct 5)'.format(57300))
ax2.set_xlabel('Days since MJD {} (2020 May 31)'.format(59000))

plt.savefig('normalised_sn_scale_change.png', bbox_inches='tight', dpi=200)
plt.savefig('normalised_sn_scale_change.eps', bbox_inches='tight')