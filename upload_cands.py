#!/usr/bin/env python

import argparse
from vcstools.client import upload_beam, upload_obsid
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
    upload_obsid(args.obsid)

    with open(args.file) as f:
        pointing_list = f.readlines()
        #for p in pointing_list:
        for i, p in enumerate(progress_bar(pointing_list, "Uploading beams")):
            name = "Blind_{}_{}".format(args.obsid, p.strip())
            upload_beam(name, args.begin, args.end,
                        mwa_search_command=args.command, mwa_search_version='v3.0', id=i)

    png_files = glob.glob('*png')
    for pf in png_files:
        print(pf)
        