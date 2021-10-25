#!/usr/bin/env python

import argparse
import glob
import concurrent.futures
import sys

from vcstools.client import upload_beam, upload_obsid, upload_cand
from vcstools.progress_bar import progress_bar
from vcstools.pointing_utils import sex2deg


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Upload the candidates on prometheus""")
    parser.add_argument('-o', '--obsid', type=str,
            help='The observation ID of the MWA observation')
    parser.add_argument('-b', '--begin', type=int,
            help="Begin time gps")
    parser.add_argument('-e', '--end', type=int,
            help="End time gps")
    parser.add_argument('-f', '--file', type=str,
            help="Pointing grid file")
    parser.add_argument('-c', '--command', type=str,
            help="mwa_search_pipeline.nf command")
    parser.add_argument('-m', '--max_threads', type=int, default=8,
            help="Maximum number of threads to use for concurrent submission.")
    parser.add_argument('-s', '--super_computer_id', type=int, default=1,
            help="Super computer ID. Ozstar: 1. Garrawarla: 2. SHAO: 3. Galaxy: 4. Magnus: 5. Default 1.")
    args = parser.parse_args()

    print("Uploading obsid")
    upload_obsid(args.obsid)

    if args.file:
        with open(args.file) as f:
            pointing_list = f.readlines()
            #for p in pointing_list:
            if args.max_threads == 1:
                # do linear method for debugging
                for i, p in enumerate(progress_bar(pointing_list, "Uploading beams")):
                    name = "Blind_{}_{}".format(args.obsid, p.strip())
                    upload_beam(name, args.begin, args.end,
                                mwa_search_command=args.command, mwa_search_version='v3.0', supercomputer_id=args.super_computer_id)
            else:
                with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_threads) as executor:
                    #for i, p in enumerate(progress_bar(pointing_list, "Uploading beams")):
                    print("Uploading beams")
                    for i, p in enumerate(pointing_list):
                        name = "Blind_{}_{}".format(args.obsid, p.strip())
                        executor.submit(upload_beam, name, args.begin, args.end,
                                        mwa_search_command=args.command, mwa_search_version='v3.0', supercomputer_id=args.super_computer_id)
    else:
        print("WARNING: Not uploading beams")

    png_files = glob.glob('B*png')
    if args.max_threads == 1:
        # do linear method for debugging
        for i, pf in enumerate(progress_bar(png_files, "Uploading {} cands".format(len(png_files)))):
            pointing = "{}_{}".format(pf.split("_")[2], pf.split("_")[3])
            # Get info from input name
            raj = pf.split("_")[2]
            decj = pf.split("_")[3]
            rad, decd = sex2deg(raj, decj)
            #print(pointing_ids[pointing], 'SMART_10min', pf, pf[:-4])
            upload_cand(rad=rad, decd=decd, obsid=args.obsid,
                        search_type='SMART_10min', search_params_id=1, png_file=pf, pfd_file=pf[:-4])
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_threads) as executor:
            print("Uploading {} cands".format(len(png_files)))
            for i, pf in enumerate(png_files):
                pointing = "{}_{}".format(pf.split("_")[2], pf.split("_")[3])
                # Get info from input name
                raj = pf.split("_")[2]
                decj = pf.split("_")[3]
                rad, decd = sex2deg(raj, decj)
                executor.submit(upload_cand, rad=rad, decd=decd, obsid=args.obsid,
                                search_type='SMART_10min', search_params_id=1, png_file=pf, pfd_file=pf[:-4])
        
