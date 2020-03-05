#!/usr/bin/env python

import argparse
from mwa_metadb_utils import get_channels


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Calculates dispersion delay""")
    parser.add_argument('-o', '--observation', type=str,
            help='The observation ID of the MWA observation')
    parser.add_argument('-lf', '--low_freq', type=float,
            help="Lowest frequency in MHz")
    parser.add_argument('-hf', '--high_freq', type=float,
            help="Highest frequency in MHz")
    parser.add_argument('-dm', type=float,
            help='DM')
    args = parser.parse_args()

    if args.observation:
        channels = get_channels(args.observation)
        lf = channels [0] * 1.28
        hf = channels [-1] * 1.28
    else:
        lf = args.low_freq
        hf = args.high_freq

    delay = 4.15*10**6 * (lf**(-2) - hf**(-2)) * args.dm
    print("Delay = {} ms".format(delay))
