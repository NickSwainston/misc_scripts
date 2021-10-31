import argparse
import numpy as np
import matplotlib.pyplot as plt


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
    import config
    from mwa_search_pipeline import process_vcs_wrapper, search_options_class
    from mwa_metadb_utils import get_channels
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
    from mwa_metadb_utils import get_channels
    import config
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
            if b'cp' in out:
                continue
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

def plot_benchmark(pns, orig, mpb_raw, orig_std, mpb_std_raw, colour, label, offset, marker,
                   read, cal, beam, write):
    # Account for the number of beams for each MPB benchmark
    mpb = []
    mpb_std = []
    for pn in pns:
        mpb.append(mpb_raw[pn-1]/pn)
        mpb_std.append(mpb_std_raw[pn-1]/pn)
    mpb = np.array(mpb)
    mpb_std = np.array(mpb_std)
    factor_improved = orig/mpb

    # Calc std
    std = []
    for i in range(len(orig_std)):
        std.append( factor_improved[i] * np.sqrt( (orig_std[i]/orig[i])**2 + (mpb_std[i]/mpb[i])**2 ) )
    
    # Plot error bars
    markersize = 3
    makerwidth = 1
    capsize = 3
    #(_, caps, _) = plt.errorbar(np.array(pns)+offset, factor_improved, yerr=std,
    #             color=colour, label=label,
    #             fmt="o", markersize=markersize, capsize=capsize)
    #for cap in caps:
    #    cap.set_markeredgewidth(makerwidth)
    plt.scatter(np.array(pns)+offset, factor_improved,
                 color=colour, label=label, marker=marker)


    # Calc theoretical improvement
    improvs = []
    for pn in pns:
        improvs.append(pn * (read + cal + beam + write) / \
                       (read + cal + pn * ( beam + write) ))
    plt.plot(np.array(pns)+offset, improvs, color=colour)

    print("Theoertical improvement for {}: {:5.2f} at {} pointings".format(label, improvs[-1], pns[-1]))
    print("Measured    improvement for {}: {:5.2f} at {} pointings".format(label, factor_improved[-1], pns[-1]))


def plot_benchmarks(max_pointings):
    pns = list(range(1,max_pointings+1))
    fig, ax = plt.subplots()

    #Ozstar MPB benchmarks serial, cal once upgrade
    ozstar_mpb_times = np.array([998.8585620833334, 692.6444897083333, 800.4549388333331, 881.5249498333333, 1006.315353, 944.5507646249999, 1065.683886, 1063.0515686249998, 1099.2300005, 1360.5995952916667, 1381.2553838333333, 1337.1792244166666, 1460.8898765000001, 1501.5252142083334, 1535.600801041667, 1561.9304510000002, 1674.5120721249998, 1750.5902402500003, 1824.18594275, 1878.5291691250002])
    ozstar_mpb_t_std = np.array([147.98231303212285, 147.24272859058829, 197.83727364838705, 171.18852669311877, 167.60135555888132, 190.65141587960483, 177.04170260993504, 169.2816422432124, 185.18820636979697, 201.6599944420951, 201.72517694475567, 138.64690917466743, 179.91234870781938, 148.1496567797671, 156.73456736977104, 185.03415973741747, 195.90517716215902, 169.5883673103956, 174.1339554201066, 120.76407871800876])
    ozstar_orig_times = np.array([972.8001408958335]*20)
    ozstar_orig_t_std = np.array([140.90075911087607]*20)

    plot_benchmark(pns, ozstar_orig_times, ozstar_mpb_times, ozstar_orig_t_std, ozstar_mpb_t_std, 'green', 'OzSTAR super computer', -0.0, ".",
                   888.9, 114.6,  42.5,  43.4)
    
    #Galaxy MPB benchmarks serial, cal once upgrade
    galaxy_mpb_times = np.array([24.946077569999996, 16.736973759999998, 13.74829847,
                        12.11295253, 11.083652268, 10.641167811666666,
                        10.711637881428569, 9.822233005000001, 9.584228442222221,
                        9.22881431, 9.166687645454546, 9.4694958725,
                        9.33200447923077, 8.810422717142856, 8.674555905999998])
    galaxy_mpb_t_std = np.array([1.4452587252036433, 0.8276700183039727, 0.7197822452137351,
                        0.38995991243788464, 0.1323128179640961, 0.38698380790094133,
                        0.23876854414689958, 0.2991264454866008, 0.28994737767526846,
                        0.2038852330143943, 0.17738044480151316, 0.20314511164648275,
                        0.19631739333354165, 0.47328326372738955, 0.19456733116853275])
    galaxy_orig_times = np.array([1.23*24]*15)
    galaxy_orig_t_std = np.array([0.03*24]*15)

    #plot_benchmark(pns, galaxy_orig_times, galaxy_mpb_times, galaxy_orig_t_std, galaxy_mpb_t_std, 'blue', 'Galaxy super computer',
    #               1.103, 0.209, 0.256, 0.062)

    #Shangia ARM MPB benchmarks serial, cal once upgrade
    #sugon gpu
    arm_mpb_times = np.array([1064.6842546666667, 1035.925497375, 1095.6567987916667, 1117.217249625, 1116.6289995, 1269.5149322916666, 1355.7624843333333, 1405.5853571666667, 1178.3422542916667, 1305.0510973333332, 1575.5143663333336, 1587.468552541667, 1656.3445657083337, 1620.6172655416667, 1714.5886795, 1890.493141375, 1839.6132066666669, 1979.3485394583338, 2031.4548157083334, 2104.5370792083336])
    arm_mpb_t_std = np.array([390.0525198905725, 358.05167142148986, 338.0677586100545, 425.76753216410646, 356.358982145573, 367.7705284620963, 341.0496397912639, 380.66996617226533, 372.9002269086173, 320.04075969979255, 442.53719150377, 436.97744920059034, 424.7550133925048, 486.53346078555603, 580.5229953920162, 509.50770291863904, 506.0167345715911, 525.3652840644522, 539.5795842738253, 559.8052841561744])
    arm_orig_times = np.array([884.2250819375]*20)
    arm_orig_t_std = np.array([318.4891566478076]*20)

    plot_benchmark(pns, arm_orig_times, arm_mpb_times, arm_orig_t_std, arm_mpb_t_std, 'red', 'CSRC prototype', 0.0, "^",
                   1329.0, 36.7, 54.5, 32.9)

    # Garrawarla 10 min test with max 20 pointings
    #pns = list(range(1,20+1))
    garra_mpb_times = np.array([490.77602183333335, 573.6949246666667, 584.4762629166667, 654.6630625833333, 692.4883127916668, 677.400559125, 741.27693525, 770.6290484583333, 787.6888258750001, 799.8723137083333, 905.1626540416668, 900.9540809583335, 942.2014737916667, 958.6411068333333, 1008.5791074999999, 1039.8349550833334, 1064.2677795416666, 1164.5895205, 1180.4333267083332, 1253.159594625])
    garra_mpb_t_std = np.array([26.513424967711305, 93.0549250040443, 86.77483618375689, 75.75238040330706, 77.90221027793791, 82.27091117305409, 83.69835235396589, 72.91378678995667, 56.3784458363581, 107.35616794927424, 85.68601358318065, 74.25646523716804, 71.9703100963701, 99.36264292041486, 67.82971562152885, 84.76268052497689, 72.0871899948102, 96.62361598666135, 89.59023889326066, 130.83673168626572])
    garra_orig_times = np.array([479.76736960416673]*20)
    garra_orig_t_std = np.array([47.1551744501299]*20)

    plot_benchmark(pns, garra_orig_times, garra_mpb_times, garra_orig_t_std, garra_mpb_t_std, 'blue', 'Garrawarla', 0.0, "*",
                   677.1, 80.6, 33.1, 20.8)



    plt.ylabel("Factor of improved processing efficiency")
    plt.xlabel("Number of simultaneous tied-array beams")
    #plt.legend(loc='upper right', bbox_to_anchor=(0.95, 0.85))
    plt.legend(loc='upper left', bbox_to_anchor=(0.005, 0.995))
    ax.set_xticks([1,5,10,15,20])
    plt.savefig("Beamformer_benchmark.eps")
    plt.savefig("Beamformer_benchmark.png", bbox_inches='tight', dpi=1000)
    

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
    parser.add_argument('--max_pointing_num', type=int, default=20,
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
    elif args.mode == 'p':
        plot_benchmarks(args.max_pointing_num)
    else:
        print("No valid option given. Doing nothing")
