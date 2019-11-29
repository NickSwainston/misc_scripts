import argparse
import numpy as np

import config
from mwa_search_pipeline import process_vcs_wrapper, search_options_class
from mwa_metadb_utils import get_channels

def make_pointing_list(pointing_in, pointing_num):
    pointing_list_list = []
    arcsec = 0
    p_ra, p_dec = pointing_in.split("_")
    dec_deg, dec_min, dec_sec = p_dec.split(":")
    for pn in range(1,pointing_num + 1):
        temp_list = []
        for n in range(1, pn + 1):
            out_dec_sec = float(dec_sec) + arcsec
            out_dec_min = int(dec_min)
            out_dec_deg = int(dec_deg)
            if out_dec_sec > 59:
                out_dec_min = out_dec_min + out_dec_sec // 60
                out_dec_sec = out_dec_sec % 60
            if out_dec_min > 59:
                out_dec_deg = out_dec_deg + out_dec_min // 60
                out_dec_min = out_dec_min % 60
            temp_list.append("{0}_{1:03d}:{2:d}:{3:05.2f}".format(p_ra,
                             out_dec_deg, int(out_dec_min), out_dec_sec))
            arcsec += 1
        pointing_list_list.append(temp_list)
    return pointing_list_list

def send_off_benchmark_jobs(obsid, cal_obs, pointing_in, begin, end, pointing_num, args, vcstools_version):
    #creating pointing list
    pointing_list_list = make_pointing_list(pointing_in, pointing_num)
    #print(pointing_list_list)
    #print(len(pointing_list_list))
    
    channels = get_channels(obsid)

    comp_config = config.load_config_file()
    DI_dir = "{0}/{1}/cal/{2}/rts/".format(comp_config['base_product_dir'],
                                           obsid, cal_obs)
    #send off jobs
    for pointing_list in pointing_list_list:
        relaunch_script = 'mwa_search_pipeline.py -o {0} -O {1} --DI_dir {2} -b {3} -e {4} --cand_type Pulsar --channels'.format(obsid, cal_obs, DI_dir, begin, end)
        for ch in channels:
            relaunch_script = "{0} {1}".format(relaunch_script, ch)
        search_opts = search_options_class(obsid, cal_id=cal_obs,
                                  begin=begin, end=end, channels=channels,
                                  DI_dir=DI_dir, relaunch_script=relaunch_script,
                                  args=args, vcstools_ver=vcstools_version)
        print("Sent off benchmarking jobs for {} pointings".format(len(pointing_list)))
        process_vcs_wrapper(search_opts, pointing_list,
                            channels=channels)
    return


def read_beanchmark_jobs(obsid, pointing, max_pointing_num, begin, end):
    comp_config = config.load_config_file()
    batch_dir = "{0}/{1}/batch/".format(comp_config['base_product_dir'], obsid)
    channels = get_channels(obsid)
    
    #get all batch files we care about
    pointing_list_list = make_pointing_list(pointing, max_pointing_num)
    import subprocess
    
    benchmark_list = []
    for pointing_list in pointing_list_list:
        pointing_str = ",".join(pointing_list)
        temp_total_time = []
        temp_read_time  = []
        temp_calc_time  = []
        temp_write_time = []
        # search for files with those pointings
        command = "grep {0} {1}mb*batch".format(pointing_str, batch_dir)
        output = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).stdout
        for out in output:
            base_batch = out.decode().split(".batch:")[0]
            begin_out = out.decode().split(" -b ")[1][:10]
            end_out   = out.decode().split(" -e ")[1][:10]
            if begin == int(begin_out) and end == int(end_out):
                with open("{0}.out".format(base_batch), "r") as batch_file:
                    lines = batch_file.readlines()
                    for line in lines:
                        if "**FINISHED BEAMFORMING**" in line:
                            temp_total_time.append(float(line.split("]")[0][1:]))
                        elif "Total read  processing time" in line:
                            temp_read_time.append(float(line.split("time: ")[1][:-3].replace(" ", "0")))
                        elif "Total calc  processing time" in line:
                            temp_calc_time.append(float(line.split("time: ")[1][:-3].replace(" ", "0")))
                        elif "Total write processing time" in line:
                            temp_write_time.append(float(line.split("time: ")[1][:-3].replace(" ", "0")))
        benchmark_list.append([len(pointing_list), temp_total_time, temp_read_time,
                               temp_calc_time, temp_write_time])
    pns = []
    total_times = []
    total_time_std = []
    for bench in benchmark_list:
        pns.append(bench[0])
        bench_total_times = bench[1]
        total_times.append(np.mean(bench_total_times)/bench[0]/100*24)
        total_time_std.append(np.std(bench_total_times)/bench[0]/100*24)

    print(f"Times: {total_times}")
    print(f"Times std: {total_time_std}")
    
    #Galaxy MPB benchmarks serial, cal once upgrade
    galaxy_mpb_times = [24.946077569999996, 16.736973759999998, 13.74829847,
                        12.11295253, 11.083652268, 10.641167811666666,
                        10.711637881428569, 9.822233005000001, 9.584228442222221,
                        9.22881431, 9.166687645454546, 9.4694958725,
                        9.33200447923077, 8.810422717142856, 8.674555905999998]
    galaxy_mpb_t_std = [1.4452587252036433, 0.8276700183039727, 0.7197822452137351,
                        0.38995991243788464, 0.1323128179640961, 0.38698380790094133,
                        0.23876854414689958, 0.2991264454866008, 0.28994737767526846,
                        0.2038852330143943, 0.17738044480151316, 0.20314511164648275,
                        0.19631739333354165, 0.47328326372738955, 0.19456733116853275]

    import matplotlib.pyplot as plt
    fig = plt.figure()
    plt.errorbar(pns, galaxy_mpb_times, yerr=galaxy_mpb_t_std, label='Galaxy multi-pixel beamformer')
    plt.errorbar(pns, [1.23*24]*15, yerr=[0.05*24]*15, label='Galaxy original beamformer')
    plt.ylabel("Processing Time per pointing per second of data (s)")
    plt.xlabel("Number of pointings")
    plt.legend(loc='upper left')
    plt.show()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Automate benchmarking""")
    parser.add_argument('-o', '--observation', type=str, default=1221832280,
            help='The observation ID of the MWA observation')
    parser.add_argument('-O', '--cal_obs', type=int, default=1221831856,
            help="Observation ID of calibrator you want to process. Used to work out "
                 "default DI_dir directory.")
    parser.add_argument('-p', '--pointing', type=str, default="21:45:50.46_-07:50:18.48",
            help='The search_opts.pointing of the fits file to be searched. The RA and Dec '
                 'seperated by _ in the format HH:MM:SS_+DD:MM:SS.')
    parser.add_argument("-b", "--begin", type=int, default=1221832351,
             help="First GPS time to process")
    parser.add_argument("-e", "--end", type=int, default=1221832450,
             help="Last GPS time to process")
    parser.add_argument('--max_pointing_num', type=int, default=15,
            help="Max number of pointings for multipixel beamformer")
    parser.add_argument('--vcstools_version', type=str, default='master',
            help="Vcstools version")
    parser.add_argument('-m', '--mode', type=str, default='s',
            help="Script mode. s: send off jobs. r: read in outputs and do statitics")
    args = parser.parse_args()

    if args.mode == 's':
        send_off_benchmark_jobs(args.observation, args.cal_obs, args.pointing,
                                args.begin, args.end, args.max_pointing_num, args,
                                args.vcstools_version)
    elif args.mode == 'r':
        dur = args.end - args.begin + 1
        read_beanchmark_jobs(args.observation, args.pointing, args.max_pointing_num,
                             args.begin, args.end)
    else:
        print("No valid option given. Doing nothing")
