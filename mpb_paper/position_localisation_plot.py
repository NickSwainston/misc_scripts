import csv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import numpy as np
from math import sqrt, exp, log, cos, sin, radians
import glob

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
#from AegeanTools.catalogs import load_catalog

from vcstools.analyse_psf import read_psf_fits

remove_average_offset = False
save_small_plots = False
add_title = False

# get data dir
import os
data_dir = os.path.dirname(os.path.realpath(__file__))[:-10] + "/data/"

# For obsid 1275758864
#fwhm = 1.05 /60. #in degrees. Estimation from vcstools scripts
# PSF elispe values from Marcin (in degrees)
emaj   = 0.02122665
emin   = 0.01499498
eangle = 154.6901
# angle at 270 degrees (for offringer_1276625432)
#fwhm = emaj * cos(radians(270)) * cos(radians(eangle)) - emin * sin(radians(270))* sin(radians(eangle))
angle_from_major = 180 - eangle
fwhm = emaj * emin / sqrt( emaj**2 * sin(radians(angle_from_major))**2 + emin**2 * cos(radians(angle_from_major))**2 )
#fwhm = radius * sin(radians(angle_from_major))

# Read in fits data
#hdul = fits.open('/astro/mwavcs/susmita/1276619416/8192x8192_1second/IDG/psf_mean.fits')
#pixel_size = 0.004 * 60 * 60 #arc seconds
#hdul = fits.open('/astro/mwavcs/nswainston/1276619416_work/PSF/pixel_0.00094deg/psf_mean.fits')

# 8.46'' pixel size data
#hdul = fits.open('/astro/mwavcs/nswainston/1276619416_work/PSF/pixel_0.00235deg/psf_mean.fits')
#offset_psf = np.array(range(25)) * pixel_size

# natural weighted psf
rav, decv, fits_data = read_psf_fits('{}wsclean_1276619416_20200619163000_briggs-psf.fits'.format(data_dir))
print(rav.shape)
print(decv.shape)
print(fits_data.shape)
fits_data = fits_data
pixel_centre = fits_data.shape[0] // 2

psf_slice = fits_data[pixel_centre][pixel_centre:pixel_centre+25]
offset_psf = (rav[pixel_centre][pixel_centre] - rav[pixel_centre][pixel_centre:pixel_centre+25]) * 60 * 60

#print(fits_data[512][512:537])
"""
print(fits_data[pixel_centre][pixel_centre:pixel_centre+25])
print(fits_data[pixel_centre][pixel_centre+4])
print(fits_data[pixel_centre+4][pixel_centre])
psf_slice = []
offest_psf = []
for i in range(pixel_centre, pixel_centre+25):
    psf_slice.append(fits_data[i][pixel_centre])
    print(rav[i][pixel_centre], rav[pixel_centre][pixel_centre])
    offest_psf.append(rav[i][pixel_centre] - rav[pixel_centre][pixel_centre])
psf_slice = np.array(psf_slice)
offest_psf = np.array(offest_psf)
"""
print(offset_psf)


print("FWHM(deg): {:.3f}".format(fwhm))
print("FWHM(``) : {:.1f}".format(fwhm*3600.))

# Read in the positions that Susmitas image using aegean
#cat = load_catalog('odo_comp.xml')
#image_positions = SkyCoord(cat, unit="deg")

# Read in data
# plot_label, cal_source, hours_away, azel_diff, file_glob
all_calibrators_data = [('offringer_1276625432',         '3C444',  1.67, 57.99, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276625432_casa_offringa/*orig_best_predicted_sn.csv'),
                        # ('offringer_imaging_1276625432', '3C444',  1.67, 57.99, None),
                        # ('RTS_1276625432',               '3C444',  1.67, 57.99, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276625432/*orig_best_predicted_sn.csv'),
                        # ('RTS_to_offringer_1276625432',  '3C444',  1.67, 57.99, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276625432_rts2ao/*orig_best_predicted_sn.csv'),
                        # ('RTS_1276468440',               '3C444', 41.9,  57.99, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276468440/*orig_best_predicted_sn.csv'),
                        # ('RTS_1276509016',               'HydA',  30.67, 46.10, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276509016/*orig_best_predicted_sn.csv'),
                        # ('RTS_1276554856',               '3C444', 17.93, 57.99, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276554856/*orig_best_predicted_sn.csv'),
                        # ('RTS_1276641272',               '3C444',  6.07, 57.99, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276641272/*orig_best_predicted_sn.csv'),
                        #No detections in the below
                        #('RTS_1276619296',      'HerA',   0.03, 23.35, '/astro/mwavcs/pulsar_search/mpb_localisation/1276619416_1276619296/*orig_best_predicted_sn.csv')
                        ]
for plot_label, cal_source, hours_away, azel_diff, file_glob in all_calibrators_data:
    # For each calibrator read in the data, calculate the offsets and other statistics and creat the plots

    if file_glob is None:
        # Susmitas imaging results
        #pulsar, psrcat_pointing, psrcat_pdmp_SN, localisation_pointing, localisation_pdmp_SN
        data = [["J1820-0427", "18:20:52.59_-04:27:37.71", 1, "18:20:52.47_-04:27:29.02", 2],
                ["J1825-0935", "18:25:30.59_-09:35:21.20" , 1, "18:25:30.52_-09:35:16.93", 2],
                ["J1833-0338", "18:33:41.89_-03:39:04.25" , 1, "18:33:41.97_-03:38:58.58", 2],
                ["J1834-0426", "18:34:25.60_-04:26:15.83" , 1, "18:34:25.55_-04:26:10.65", 2],
                ['J1913-0440', "19:13:54.17_-04:40:47.68" , 1, "19:13:54.08_-04:40:46.25", 2]]
    else:
        data = []
        csv_file_list = glob.glob(file_glob)
        csv_file_list.sort()
        for file_loc in csv_file_list:
            with open(file_loc) as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                csv_data = []
                for row in spamreader:
                    csv_data.append(row)
                    #print(row)
                pulsar = file_loc.split('/')[-1].split('_')[2]
                #print(csv_data)
                #pulsar, psrcat_pointing, psrcat_pdmp_SN, localisation_pointing, localisation_pdmp_SN
                data.append([pulsar, csv_data[0][0], float(csv_data[0][1]), csv_data[1][0], float(csv_data[1][1])])


    # Calculate offsets and improvements
    positions_calcs = []
    pulsars = []
    distances = []
    brightest_sns = []
    for cand in data:
        pulsar, psrcat_pointing, psrcat_pdmp_SN, localisation_pointing, localisation_pdmp_SN = cand

        """
        if plot_label == 'offringer_1276625432':
            # Find image sources within positions
            offset_search_radius = 0.4 #arcminute
            print("{} tied array beamforming".format(pulsar))
            loc_point = SkyCoord( localisation_pointing.split("_")[0], localisation_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
            for ip in image_positions:
                if ip.separation(loc_point   ).arcminute < offset_search_radius:
                    print("        {:6.2f} {:6.2f} {}".format(ip.ra.deg, ip.dec.deg, ip.separation(loc_point).arcsec))
            print("{} catalogue".format(pulsar))
            psrcat_point = SkyCoord( psrcat_pointing.split("_")[0], psrcat_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
            for ip in image_positions:
                if ip.separation(psrcat_point).arcminute < offset_search_radius:
                    print("        {:6.2f} {:6.2f} {}".format(ip.ra.deg, ip.dec.deg, ip.separation(psrcat_point).arcsec))
            print("")
        """


        pulsars.append(pulsar)
        improvement = ( localisation_pdmp_SN / psrcat_pdmp_SN - 1 ) * 100
        if float(improvement) > 0.:
            c1 = SkyCoord( psrcat_pointing.split("_")[0], psrcat_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
            og_ra  = c1.ra.deg
            og_dec = c1.dec.deg

            c2 = SkyCoord( localisation_pointing.split("_")[0], localisation_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
            loc_ra  = c2.ra.deg
            loc_dec = c2.dec.deg

            ra_diff = loc_ra - og_ra
            dec_diff = loc_dec - og_dec
            distance  = c1.separation(c2).deg
            #print("diff(ra dec ``): {:6.2f} {:6.2f}".format(ra_diff*60.*60., dec_diff*60.*60.,))
            positions_calcs.append([pulsar, og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement, distance])

            # work out expected improvement
            #distance = np.sqrt(cos(radians(loc_dec))**2 * ra_diff**2 + dec_diff**2)
            distances.append(distance)
            gaussian = np.exp(-(distance)**2/((fwhm*0.85)**2))
            improvement_expected = 100 * (1 - gaussian)
            #print("")
            #print("Pulsar: {}  distance(``): {:.1f}".format(pulsar, distance*3600))
            #print("expected SN improvement(%): {:.2f}".format(improvement_expected))
            #print("calculated  improvement(%): {}".format(improvement))
            #print(og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement)
            if pulsar in ["J1820-0427", "J1823+0550", "J1825-0935", "J1833-0338", "J1834-0426", "J1913-0440"]:
                #print(pulsar, psrcat_pointing, psrcat_pdmp_SN, localisation_pointing, localisation_pdmp_SN)
                brightest_sns.append(localisation_pdmp_SN)
        else:
            c = SkyCoord( psrcat_pointing.split("_")[0], psrcat_pointing.split("_")[1], frame='icrs', unit=(u.hourangle,u.deg))
            og_ra  = c.ra.deg
            og_dec = c.dec.deg
            distances.append(0.)
            #positions_calcs.append([pulsar, og_ra, og_dec, og_ra, og_dec, 0., 0., 0.])
            if pulsar in ["J1820-0427", "J1823+0550", "J1825-0935", "J1834-0426", "J1913-0440"]:
                brightest_sns.append(psrcat_pdmp_SN)

        #if plot_label == 'offringer_1276625432':
        #    print(f"{distance},{improvement}")
    #print(np.array(positions_calcs))

    # Save data
    #np.savetxt("{}.txt".format(plot_label), np.array(positions_calcs), delimiter=',', fmt='%s',
    #           header="pulsar, psrcat RA (deg), psrcat Dec (deg), localisation RA (deg), localisation Dec (deg), RA diff (deg), Dec diff (deg), SN improvement (%), Offset (deg)")

    # work out average dec and ra diffs
    ra_off_mean = np.mean(np.array(np.array(positions_calcs)[:,5], dtype=float))
    dec_off_mean = np.mean(np.array(np.array(positions_calcs)[:,6], dtype=float))


    # Plot data
    colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k',
            'tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan',
            'indigo', 'navy', 'chocolate']
    size = 5

    # 3 plot method
    #bigfig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(3*size,size))

    # inset plot method
    bigfig, (ax2, ax3) = plt.subplots(1,2, figsize=(2*size,size))
    ax1 = plt.axes([0,0,1,1])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax2, [0.64,0.1,0.3,0.3])
    ax1.set_axes_locator(ip)

    # First plot. Relative offsets
    plt.figure(figsize=(size, size))
    for ci, pos in enumerate(positions_calcs):
        pulsar, og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement, distance = pos
        if remove_average_offset and 'imaging' not in plot_label:
            ra_diff = ra_diff - ra_off_mean
            dec_diff = dec_diff - dec_off_mean
        #if pulsar in ["J1820-0427", "J1825-0935", "J1833-0338", "J1834-0426", "J1913-0440"]:
        plt.arrow(0., 0., ra_diff*60.*60., dec_diff*60.*60., head_width=1.5, head_length=2.5, fc=colours[ci], ec=colours[ci], label=pulsar)
        ax1.arrow(0., 0., ra_diff*60.*60., dec_diff*60.*60., head_width=1.5, head_length=2.5, fc=colours[ci], ec=colours[ci], label=pulsar)
    pos_size = 60
    plt.xlim(-pos_size, pos_size)
    ax1.set_xlim(-pos_size, pos_size)
    plt.ylim(-pos_size, pos_size)
    ax1.set_ylim(-pos_size, pos_size)
    ax1.set_xticklabels([0, r"-50$^{\prime\prime}$", "0", r"50$^{\prime\prime}$"])
    ax1.set_yticklabels([0, r"-50$^{\prime\prime}$", "0", r"50$^{\prime\prime}$"])
    #plt.legend(loc='upper right')
    #plt.xlabel(r"$\Delta$ Right Ascension ($^{\prime\prime}$)")
    #ax1.set_xlabel(r"$\Delta$ Right Ascension ($^{\prime\prime}$)")
    #plt.ylabel(r"$\Delta$ Declination ($^{\prime\prime}$)")
    #ax1.set_ylabel(r"$\Delta$ Declination ($^{\prime\prime}$)")
    if save_small_plots:
        plt.savefig("{}_relative_offsets.png".format(plot_label), dpi=500)
    else:
        plt.close()

    # Plot relative positions on the sky
    plt.figure(figsize=(size, size))
    for ci, pos in enumerate(positions_calcs):
        pulsar, og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement, distance = pos
        scale = 100.
        #print(og_ra, og_dec, ra_diff*scale, dec_diff*scale)
        plt.arrow(og_ra, og_dec, ra_diff*scale, dec_diff*scale, head_width=distance*0.05, head_length=distance*0.1, fc=colours[ci], ec=colours[ci])
        ax2.arrow(og_ra, og_dec, ra_diff*scale, dec_diff*scale, head_width=distance*0.05, head_length=distance*0.1, fc=colours[ci], ec=colours[ci])
    xlim_min, xlim_max, ylim_min, ylim_max = [260., 295., -15., 2.]
    #xlim_min, xlim_max, ylim_min, ylim_max = [150., 195., -36., -53.]
    plt.xlim(xlim_min, xlim_max)
    ax2.set_xlim(xlim_min, xlim_max)
    plt.ylim(ylim_min, ylim_max)
    ax2.set_ylim(ylim_min, ylim_max)

    ax2.set_xticklabels(["15:20", "15:40", "16:00", "16:20",
                         "16:40", "17:00", "17:20", "17:40"])

    plt.xlabel(r"Right Ascension (HH:MM)")
    ax2.set_xlabel(r"Right Ascension (HH:MM)")
    plt.ylabel(r"Declination ($^\circ$)")
    ax2.set_ylabel(r"Declination ($^\circ$)")
    #plt.legend(loc='upper left')
    if save_small_plots:
        plt.savefig("{}_sky_offsets.png".format(plot_label), dpi=500)
    else:
        plt.close()

    # Plot relative positions on the sky
    plt.figure(figsize=(size, size))
    legend_lines = []
    for ci, pos in enumerate(positions_calcs):
        pulsar, og_ra, og_dec, loc_ra, loc_dec, ra_diff, dec_diff, improvement, distance = pos
        distance = distance * 3600
        #if pulsar in ["J1820-0427", "J1825-0935", "J1834-0426", "J1913-0440"]:
        plt.scatter(distance, improvement, c=colours[ci])
        legend_lines.append(ax3.scatter(distance, improvement, c=colours[ci]))
    # Add expected improvement gaussian
    x = np.linspace(0, pos_size, 120)
    fwhm_as = fwhm * 60 * 60
    sigma = fwhm_as / 2.35482
    gaussian = np.exp(-(x)**2/(2*sigma**2))
    #print(gaussian)
    improvement_expected = 100 * (1 - gaussian)
    #plt.plot(x, improvement_expected)
    #legend_lines.append(ax3.plot(x, improvement_expected))

    plt.plot(offset_psf, fits_data[pixel_centre][pixel_centre:pixel_centre+25])
    #legend_lines.append(ax3.plot(offset_psf, 100 - fits_data[pixel_centre][pixel_centre:pixel_centre+25]*100))
    legend_lines.append(ax3.plot(offset_psf, 100 - psf_slice*100))

    plt.xlabel(r"Offset ($^{\prime\prime}$)")
    ax3.set_xlabel(r"Offset ($^{\prime\prime}$)")
    plt.ylabel("SN Degredation (%)")
    ax3.set_ylabel("SN Degredation (%)")
    plt.xlim(0, pos_size)
    ax3.set_xlim(0, pos_size)
    plt.ylim(60, 0)
    ax3.set_ylim(60, 0)
    #plt.legend(loc='upper left')
    if save_small_plots:
        plt.savefig("{}_SN_improvements_vs_offsets.png".format(plot_label), dpi=500)
    else:
        plt.close()

    # Combined plot
    bigfig.legend(
        handles=legend_lines,     # The line objects
        labels=["Exp improve"] + pulsars,#["J1820-0427", "J1825-0935", "J1834-0426", "J1913-0440"],#pulsars,   # The labels for each line
        loc="center right",   # Position of legend
        borderaxespad=0.1,    # Small spacing around legend box
        title="Pulsars",  # Title for the legend
        framealpha=1.0,
    )
    if add_title:
        bigfig.suptitle(f"{plot_label}  {cal_source}   mean offset: {np.mean(distances)*60*60:.2f}''   mean brightest SN {np.mean(brightest_sns):.2f}   Hrs Away: {hours_away}   Az El diff: {azel_diff}$^\circ$")
    #bigfig.legend(handles, labels, loc='upper left')
    #bigfig.subplots_adjust(right=0.85)
    #bigfig.axis('equal')
    bigfig.savefig("{}_all_offsets.pdf".format(plot_label), format='pdf')
