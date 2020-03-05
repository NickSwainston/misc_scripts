import argparse
import subprocess
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
import re

import os
import psrqpy
import pandas

import prof_utils
from prof_utils import est_sn_from_prof

import matplotlib.pyplot as plt

def plot_sn(obsid, pulsar):
    if pulsar.startswith("J"):
        pulsar = pulsar[1:]
    bestprof_files = []
    pointings = os.listdir("/group/mwaops/vcs/{0}/pointings".format(obsid)) 
    #print(pointings)
    for point in pointings:
        #file_name = "/group/mwaops/vcs/{0}/pointings/{1}/{1}_PSR_{2}.pfd.bestprof".format(obsid, point, pulsar)
        file_name = "/group/mwaops/vcs/{0}/pointings/{1}/{0}_{1}_PSR_{2}.pfd.bestprof".format(obsid, point, pulsar)
        #print(file_name)
        if os.path.exists(file_name):
            bestprof_files.append([file_name, point])
    #print(raw_grep)
    
    data = []
    max_sn = 0.
    for d in bestprof_files:
        bestprof_data = prof_utils.get_from_bestprof(d[0])
        obsid, pulsar, dm, period, period_uncer, obsstart, time_detection,\
                                       profile, num_bins = bestprof_data
        print(d)
        try:
            sn = est_sn_from_prof(profile, float(period))[0]
        except:
            sn = 0
        if sn is None:
            sn = 0
        if max_sn < sn:
            max_sn = sn
        data.append([d[1], sn])
    
    #get pulsar position
    if not pulsar.startswith("J"):
        pulsar = "J" + pulsar
    query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=os.environ['PSRCAT_FILE']).pandas
    raj = query["RAJ"][0]
    decj = query["DECJ"][0]
    
    fig = plt.figure()
    ras = []
    decs = []
    coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
    cat_ra = coord.ra.degree
    cat_dec = coord.dec.degree
    for d in data:
        pointing, sn = d
        coord = SkyCoord(pointing.split("_")[0], pointing.split("_")[1],unit=(u.hourangle,u.deg))
        ra = coord.ra.degree - cat_ra
        dec = coord.dec.degree - cat_dec
        #print(ra, dec, sn)
        ras.append(ra)
        decs.append(dec)
        if sn > max_sn / 2.:
            colour = 'orange'
        else:
            colour = 'black'
        if sn == max_sn:
            colour = 'red'
        plt.text(ra, dec, "{:.2f}".format(sn), va='center', ha='center', fontsize=6, color=colour)
    print([min(ras)-0.01, max(ras)+0.01, min(decs)-0.01, max(decs)+0.01])
    plt.axis([min(ras)-0.01, max(ras)+0.01, min(decs)-0.01, max(decs)+0.01])
    plt.xlabel(r"$\Delta$ Right Ascension ($^{\circ}$)")
    plt.ylabel(r"$\Delta$ Declination ($^{\circ}$)")
    plt.savefig("{0}_{1}_sn.png".format(obsid, pulsar), dpi=1000)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
             For a cross pattern of pointings (probably done by grid.py -m cross)
             that have been prepfolded, uses their signal to noise to work out the fwhm.
             """)
    parser.add_argument('-o','--obsid',type=str,help='The observation ID to be searched')
    parser.add_argument('-p','--pulsar',type=str,help='The pulsar.')
    args=parser.parse_args()
    plot_sn(args.obsid, args.pulsar)
