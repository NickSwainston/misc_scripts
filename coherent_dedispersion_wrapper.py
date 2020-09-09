import psrqpy
import pandas
import numpy as np

from job_submit import submit_slurm

def dms_47_tuc():
    params = ['JNAME', 'DM']
    pulsar_list = ["J0024-7204aa", "J0024-7204ab", "J0024-7204C", "J0024-7204D",
                   "J0024-7204E", "J0024-7204F", "J0024-7204G", "J0024-7204H",
                   "J0024-7204I", "J0024-7204J", "J0024-7204L", "J0024-7204M",
                   "J0024-7204N", "J0024-7204O", "J0024-7204P", "J0024-7204Q",
                   "J0024-7204R", "J0024-7204S", "J0024-7204T", "J0024-7204U",
                   "J0024-7204V", "J0024-7204W", "J0024-7204X", "J0024-7204Y",
                   "J0024-7204Z"]
    query = psrqpy.QueryATNF(params=params, psrs=pulsar_list).pandas

    for pi, _ in enumerate(pulsar_list):
        pulsar = query["JNAME"][pi]
        dm = query["DM"][pi]
        if pulsar == "J0024-7204V":
            dm_err = 0.008
        else:
            dm_err = query["DM_ERR"][pi]
        print(pulsar,dm,dm_err)
        #for dm_step in np.arange(dm-dm_err, dm+dm_err, 0.000014744):
        #for dm_step in np.arange(dm-dm_err, dm+dm_err, 0.00005):
        for dm_step in [dm]:

            #send off job
            commands = []
            commands.append("cd /group/mwaops/vcs/1227009976/pointings/00:24:07.01_-72:04:32.91")
            commands.append("psrcat -e {0} > {0}.par".format(pulsar))
            commands.append("mkdir -p {0}".format(pulsar))
            commands.append("cd {0}".format(pulsar))
            commands.append("for CHAN in $(seq 109 132); do")
            commands.append("   #dspsr -b 100 -O {0}_DM_{1}_ch${{CHAN}} -E {0}.par -A -L 10 -cont -U 590 G0057_1227009976_ch${{CHAN}}_u.hdr".format(pulsar, dm_step))
            commands.append("   mkdir -p ch_$CHAN")
            commands.append("   cd ch_$CHAN")
            commands.append("   ln -s /group/mwaops/vcs/1227009976/pointings/00:24:07.01_-72:04:32.91/G0057_1227009976_ch${CHAN}_u.hdr G0057_1227009976_ch${CHAN}_u.hdr")
            commands.append("   ln -s /group/mwaops/vcs/1227009976/pointings/00:24:07.01_-72:04:32.91/G0057_1227009976_ch${CHAN}_u.vdif G0057_1227009976_ch${CHAN}_u.vdif")
            commands.append("   dspsr -b 128 -E ../../{0}.par -cont -U 590 -J ../../giant_pulse.psh -s -K G0057_1227009976_ch${{CHAN}}_u.hdr".format(pulsar, dm_step))
            commands.append("   cd ..")
            commands.append("done")
            commands.append("#psradd -R -m time {0}_DM_{1}_ch*.ar -o {0}_DM_{1}.ar".format(pulsar, dm_step))
            commands.append("#pav -DTFCp -g {0}_DM_{1}.prof.ps/cps {0}_DM_{1}.ar".format(pulsar, dm_step))
            submit_slurm("{0}_DM_{1}".format(pulsar, dm_step), commands,
                 batch_dir="/group/mwaops/vcs/1227009976/batch/tuc_search",
                 slurm_kwargs={"time": "23:00:00"},
                 module_list=['dspsr'],
                 submit=True, mem="1024")

if __name__ == '__main__':
    dms_47_tuc()
