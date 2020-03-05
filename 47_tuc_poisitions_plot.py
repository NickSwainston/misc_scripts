import argparse
import subprocess
import os
import psrqpy
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import pandas
import statistics 
import math

def plot_positions():
    pulsar = ["J0024-7204aa", "J0024-7204ab", "J0024-7204C", "J0024-7204D", "J0024-7204E", "J0024-7204F", "J0024-7204G", "J0024-7204H", "J0024-7204I", "J0024-7204J", "J0024-7204L", "J0024-7204M", "J0024-7204N", "J0024-7204O", "J0024-7204P", "J0024-7204Q", "J0024-7204R", "J0024-7204S", "J0024-7204T", "J0024-7204U", "J0024-7204V", "J0024-7204W", "J0024-7204X", "J0024-7204Y", "J0024-7204Z"]
    query = psrqpy.QueryATNF(psrs=pulsar).pandas 
    ras = []
    decs = []
    ax = plt.gca()
    for pi, _ in enumerate(pulsar):
        raj = query["RAJ"][pi]
        decj = query["DECJ"][pi]
        coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
        ras.append(coord.ra.degree)
        decs.append(coord.dec.degree)
    plt.scatter(ras, decs)

    #work out where a circle that contains all pulsars
    avg_ra  = statistics.mean(ras)
    avg_dec = statistics.mean(decs)
    max_dis = 0.
    for pi, _ in enumerate(pulsar):
        dis = math.sqrt( (ras[pi]-avg_ra)**2 + (decs[pi]-avg_dec)**2)
        if dis > max_dis:
            max_dis = dis

    plt.scatter(avg_ra, avg_dec, color='r')
    circle1 = plt.Circle((avg_ra, avg_dec), max_dis , color='r', fill=False)
    ax.add_artist(circle1)

    plt.xlabel(r"Right Ascension ($^{\circ}$)")
    plt.ylabel(r"Declination ($^{\circ}$)")
    plt.savefig("47_tuc_positions.eps")

    c = SkyCoord( avg_ra, avg_dec, frame='icrs', unit=(u.deg,u.deg))
    rajs = c.ra.to_string(unit=u.hour, sep=':')
    decjs = c.dec.to_string(unit=u.degree, sep=':')
    print("Centre: {0}_{1}".format(rajs, decjs))
    print("Radius: {0} deg".format(max_dis))

if __name__ == '__main__':
    plot_positions()
