#!/usr/bin/env python

import argparse
import glob
from vcstools.client import upload_beam, upload_obsid, upload_cand
from vcstools.progress_bar import progress_bar


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
    args = parser.parse_args()

    print("Uploading obsid")
    #upload_obsid(args.obsid)

    pointing_ids = {}
    with open(args.file) as f:
        pointing_list = f.readlines()
        #for p in pointing_list:
        for i, p in enumerate(progress_bar(pointing_list, "Uploading beams")):
            name = "Blind_{}_{}".format(args.obsid, p.strip())
            pointing_ids[p.strip()] = i + 4
            #upload_beam(name, args.begin, args.end,
            #            mwa_search_command=args.command, mwa_search_version='v3.0', supercomputer_id=1)
    png_files = glob.glob('B*png')
    for i, pf in enumerate(progress_bar(png_files, "Uploading cands")):
        pointing = "{}_{}".format(pf.split("_")[2], pf.split("_")[3])
        #print(pointing_ids[pointing], 'SMART_10min', pf, pf[:-4])
        upload_cand(pointing_ids[pointing], search_type='SMART_10min', search_params_id=1, png_file=pf, pfd_file=pf[:-4])
        
