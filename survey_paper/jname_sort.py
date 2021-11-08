import csv
import os
import glob
import logging
import psrqpy
import subprocess
from shutil import copyfile

from vcstools.pointing_utils import format_ra_dec

logger = logging.getLogger(__name__)

def get_db_auth_addr():
    """
    Checks for MWA database usernames and passwords

    Returns:
    --------
    auth: tuple
        The username and password for the pulsar databse
    web_address: string
        The web address of the pulsar database
    """
    web_address = 'https://pulsar-cat.icrar.uwa.edu.au/'
    if 'MWA_PULSAR_DB_USER' in os.environ and 'MWA_PULSAR_DB_PASS' in os.environ:
        auth = (os.environ['MWA_PULSAR_DB_USER'],os.environ['MWA_PULSAR_DB_PASS'])
    else:
        auth = None
        raise NoAuthError(  """
                            No MWA Pulsar Database username/password found
                            Please add the following to your .bashrc:
                            'export MWA_PULSAR_DB_USER="<username>"'
                            'export MWA_PULSAR_DB_PASS="<password>"'
                            'replacing <username> <password> with your MWA Pulsar Database username and password
                            """ )

    return web_address, auth

submit = True
calid_dict = {'1302282040': 1302171592, '1302540536': 1302603296, '1302712864': 1302732080, '1302106648': 1302085256, '1301847296': 1301826240, '1301674968': 1301739904, '1301412552': 1301480888, '1301240224': 1301221888, '1300981728': 1300962880, '1300809400': 1300790208, '1268321832': 1268321544, '1268063336': 1268063056, '1267459328': 1267459048, '1267283936': 1267268464, '1267111608': 1267111328, '1266680784': 1266664072, '1265470568': 1265470280, '1264867416': 1264850632, '1266932744': 1266932464, '1266329600': 1266329312, '1266155952': 1266155672, '1265725128': 1265724848, '1265983624': 1265983344, '1261241272': 1261240992, '1260638120': 1260637840, '1259427304': 1259427008, '1259685792': 1259691192, '1258221008': 1258219696, '1257617424': 1257622816, '1257010784': 1257010488, '1256407632': 1256381240, '1255803168': 1255899632, '1255197408': 1255208880, '1254594264': 1254593968, '1253991112': 1254000168, '1253471952': 1253477344, '1252780888': 1252705176, '1252177744': 1252100856, '1225462936': 1225478264, '1225118240': 1225133000, '1224859816': 1224874512, '1227009976': 1227007528, '1224252736': 1224277624, '1255444104': 1255443816, '1226062160': 1226054696, '1225713560': 1225710296, '1223042480': 1223068840, '1222697776': 1222695592, '1222435400': 1222434520, '1221832280': 1221831856, '1221399680': 1221342176}

# grabed from the R column of the SMART google sheet
jname_str = "J2241-5236 J2145-0750 J2222-0137 J2248-0101 J2155-3118 J2048-1616 J2108-3429 J2330-2005 J2325-0530 J0034-0721 J2145-0750 J2336-01 J2234+2114 J2317+2149 J2241-5236 J2330-2005 J0152-1637 J0206-4028 J2354-22 J0134-2937 J0038-2501 J0030+0451 J0034-0721 J0034-0534 J2241-5236 J0133-6957 J2324-6054 J0206-4028 J0255-5304 J0051+0423 J0152-1637 J0151-0635 J0152-1637 J0255-5304 J0206-4028 J0418-4154 J0304+1932 J0450-1248 J0452-1759 J0459-0210 J0401-7608 J0450-1248 J0459-0210 J0452-1759 J0520-2553 J0450-1248 J0514-4408 J0636-4549 J0702-4956 J0600-5756 J0437-4715 J0630-2834 J0737-3039A J0636-4549 J0534+2200 J0528+2200 J0450-1248 J0452-1759 J0459-0210 J0601-0527 J0624-0424 J0630-2834 J0820-1350 J0737-3039A J0742-2822 J0729-1836 J0742-2822 J0749-4247 J0820-4114 J0820-3921 J0835-4510 J0837-4135 J0838-3947 J0924-5302 J0942-5552 J0955-5304 J0959-4809 J0826+2637 J0837+0610 J0922+0638 J0729-1448 J0729-1836 J0742-2822 J0758-1528 J0820-1350 J0930-2301 J0837-4135 J0856-6137 J0905-6019 J0924-5302 J0924-5814 J0942-5657 J0842-4851 J0907-5157 J0902-6325 J0856-6137 J0902-6325 J0904-7459 J0905-6019 J0924-5302 J0924-5814 J0942-5552 J0942-5657 J0922+0638 J1022+1001 J0953+0755 J0908-1739 J0820-4114 J0742-2822 J0820-1350 J0835-4510 J0837-4135 J0855-3331 J1012-2337 J0835-4510 J0837-4135 J0924-5302 J0942-5657 J0955-5304 J0959-4809 J1003-4747 J1057-5226 J1116-4122 J0922+0638 J0943+1631 J0946+0951 J0953+0755 J0908-1739 J0944-1354 J1018-1642 J1057-5226 J1059-5742 J1116-4122 J1121-5444 J1123-4844 J1123-6651 J1136-5525 J1141-6545 J1146-6030 J1202-5820 J1224-6407 J1225-5556 J1240-4124 J1312-5402 J1320-5359 J0953+0755 J1136+1551 J1012-2337 J1018-1642 J1034-3224 J1041-1942 J1311-1228 J1257-1027 J1313+0931 J1300+1240 J1059-5742 J1202-5820 J1224-6407 J1430-6623 J1141-6545 J1112-6926 J1239-6832 J1456-6843 J1340-6456 J1123-6651  J1328-4357 J1240-4124 J1355-5153 J1418-3921 J1335-3642 J1320-5359 J1311-1228 J1332-3032 J1418-3921 J1313+0931 J1300+1240 J1311-1228 J1453-6413 J1456-6843 J1430-6623 J1534-5334 J1440-6344 J1355-5153 J1418-3921 J1543-0620 J1543+0929 J1510-4422 J1527-3931 J1455-3330 J1507-4352 J1418-3921 J1534-5334 J1536-4948"
jname_str = jname_str.replace("  ", " ")
jname_list = list(dict.fromkeys(jname_str.split(" ")))
jname_dict = {}
for jn in jname_list:
    jname_dict[jn] = []

with open("{}/smart.csv".format(os.path.dirname(os.path.realpath(__file__))), "r") as it:
    reader = csv.reader(it)
    csv_data = []
    for row in reader:
        csv_data.append(row)

for o, o_jnames in csv_data:
    o_jnames = o_jnames.split(" ")
    for oj in o_jnames:
        jname_dict[oj].append(o)

obs_with_missing_profiles = []
query = psrqpy.QueryATNF(params = ['PSRJ', 'DM', 'P0', 'RAJ', 'DECJ']).pandas
with open('SMART_pulsars.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile)
    for jname, obsids in sorted(jname_dict.items()):
        for o, obsid in enumerate(obsids):
            if o == 0:
                joutput = jname
            else:
                joutput = ""
            query_id = list(query['PSRJ']).index(jname)
            # See if dpp ran and what is the best bin size
            #print("/astro/mwavcs/vcs/{0}/dpp/{0}_{1}/{0}_{1}*bins_gaussian_fit.png".format(
            #                          obsid, jname))

            # Finding on Garra and pulsar database
            """
            gaussian_fits = glob.glob("/astro/mwavcs/vcs/{0}/dpp/{0}_{1}/{0}_{1}*bins_gaussian_fit.png".format(
                                      obsid, jname))
            bestprof = None
            png = None
            archive = None

            if len(gaussian_fits) != 0:
                # Grab best bin profile
                all_bins = []
                for g in gaussian_fits:
                    all_bins.append(int(g.split("_bins_gaussian_fit.png")[0].split("_")[-1]))
            else:
                # see what bin profiles there are
                dpp_bestprof = glob.glob("/astro/mwavcs/vcs/{0}/dpp/{0}_{1}/pf_{0}_{1}_*_b*bestprof".format(
                                         obsid, jname))
                all_bins = []
                for b in dpp_bestprof:
                    all_bins.append(int(b.split("_b")[-1].split("_")[0]))
            if all_bins:
                bins = max(all_bins)
                if bins < 100:
                    # MSP so use smallest bin profile
                    bins = min(all_bins)
                
                # Find bestprof
                #print("/astro/mwavcs/vcs/{0}/dpp/{0}_{1}/pf_{0}_{1}_*_b{2}*bestprof".format(
                #                         obsid, jname, bins))
                dpp_bestprof = glob.glob("/astro/mwavcs/vcs/{0}/dpp/{0}_{1}/pf_{0}_{1}_*_b{2}*bestprof".format(
                                        obsid, jname, int(bins)))
                if len(dpp_bestprof) != 0:
                    bestprof = dpp_bestprof[0]

                # Find png
                dpp_png = glob.glob("/astro/mwavcs/vcs/{0}/dpp/{0}_{1}/pf_{0}_{1}_*_b{2}*png".format(
                                        obsid, jname, int(bins)))
                if len(dpp_png) != 0:
                    png = dpp_png[0]
            
                # Find archive 
                dpp_archive = glob.glob("/astro/mwavcs/vcs/{0}/dpp/{0}_{1}/{0}_{1}_archive.ar".format(
                                        obsid, jname))
                if len(dpp_archive) == 0:
                    print("No archive for {} {}".format(obsid, jname))
                else:
                    archive = dpp_archive[0]
                    print(archive)

            if not png or not bestprof:
                # Attempt to download the data from the pulsar database
                from mwa_pulsar_client.client import detection_file_download
                web_address, auth = get_db_auth_addr()

                if not os.path.isdir("/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{}".format(obsid)):
                    os.mkdir("/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{}".format(obsid))
                if not png:
                    # download ps
                    possible_names = ["{}_{}.prepfold.ps".format(obsid, jname),
                                      "{}_{}.prepfold.ps".format(obsid, jname.replace("+","")),
                                      "{}_{}_c1221342176_b22.prepfold.ps".format(obsid, jname),
                                      "{}_{}_c1268063056_b1024.prepfold.ps".format(obsid, jname.replace("+",""))]
                    for pn in possible_names:
                        if os.path.isfile("/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{0}/{1}".format(obsid, pn)):
                            png = "/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{0}/{1}".format(obsid, pn)
                            break
                        try:
                            detection_file_download(web_address, auth,
                                                    pn,
                                                    "/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{}".format(obsid))
                            png = "/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{0}/{0}_{1}.prepfold.ps".format(obsid, jname)
                        except:
                            logger.debug("{} not on database".format(pn))
                
                if not bestprof:
                    # download bestprof
                    possible_names = ["{}_{}.bestprof".format(obsid, jname),
                                      "{}_{}.bestprof".format(obsid, jname.replace("+","")),
                                      "{}_{}_c1221342176_b22.bestprof".format(obsid, jname),
                                      "{}_{}_c1268063056_b1024.bestprof".format(obsid, jname.replace("+",""))]
                    for pn in possible_names:
                        if os.path.isfile("/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{0}/{1}".format(obsid, pn)):
                            bestprof = "/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{0}/{1}".format(obsid, pn)
                            break
                        try:
                            detection_file_download(web_address, auth,
                                                    pn,
                                                    "/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{}".format(obsid))
                            bestprof = "/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{0}/{0}_{1}.bestprof".format(obsid, jname)
                        except:
                            logger.debug("{} not on database".format(pn))
            """

            # Find in SMART quick look directory (must be run from it)
            bestprof = glob.glob("{0}/*{0}_{1}*f".format(obsid, jname))
            if len(bestprof) == 1:
                bestprof = bestprof[0]
            png = glob.glob("{0}/*{0}_{1}*ps".format(obsid, jname)) + glob.glob("{0}/*{0}_{1}*png".format(obsid, jname))
            if len(png) == 1:
                png = png[0]
            sefd_loc = glob.glob("{0}/{1}_{0}*stats".format(obsid, jname))
            # Grab DM and SN
            if bestprof:
                with open(bestprof, "r") as bestprof_f:
                    lines = bestprof_f.readlines()
                    dm = lines[14][22:-1]
                    sigma = lines[13].split("(~")[-1].split(" sigma")[0]
            else:
                dm = None
                sigma = None

            # output results
            spamwriter.writerow([joutput, obsid, query["P0"][query_id], query["DM"][query_id], dm, sigma, bestprof, png])
            if not png or not bestprof:
                print("{},{},{},{}".format(joutput, obsid, bestprof, png))
                obs_with_missing_profiles.append(obsid)
            elif submit:
                ra_dec_list = format_ra_dec([[query["RAJ"][query_id], query["DECJ"][query_id]]])
                pointing = "{}_{}".format(ra_dec_list[0][0], ra_dec_list[0][1])
                #print("/astro/mwavcs/vcs/{0}/sefd_simulations/{1}_{0}_*.stats".format(obsid, jname))
                #sefd_loc = glob.glob("/astro/mwavcs/vcs/{0}/sefd_simulations/{1}_{0}_*.stats".format(obsid, jname))
                if len(sefd_loc) == 0:
                    print("Missing sefd for {} {} so resubmitting".format(obsid, jname))
                    command = "submit_to_database.py -o {} -O {} -b {} -p {} --pointing {} --vcstools_version nswainston --dont_upload".format(obsid, calid_dict[obsid], bestprof, jname, pointing)
                    #print(command)
                    #test = subprocess.Popen(command.split(" "), stdout=subprocess.PIPE)
                    #output = test.communicate()[0]
                    #os.system(command)
                else:
                    print("Found sefd for {} {} so running flux calc".format(obsid, jname))
                    command = "submit_to_database.py -o {} -O {} -b {} -p {} --pointing {} --vcstools_version nswainston --dont_upload --sefd_file {}".format(obsid, calid_dict[obsid], bestprof, jname, pointing, sefd_loc[0])
                    print(command)
                    os.system(command)


                
                # move files to saved dir
                """
                dst = "/astro/mwavcs/pulsar_search/SMART_quick_look_detection/{}".format(obsid)
                if not os.path.isdir(dst):
                    os.mkdir(dst)
                if not os.path.isfile("{}/{}".format(dst, png.split("/")[-1])):
                    copyfile(png, "{}/{}".format(dst, png.split("/")[-1]))
                if not os.path.isfile("{}/{}".format(dst, bestprof.split("/")[-1])):
                    copyfile(bestprof, "{}/{}".format(dst, bestprof.split("/")[-1]))
                if not os.path.isfile("{}/{}".format(dst, sefd_loc[0].split("/")[-1])):
                    copyfile(sefd_loc[0], "{}/{}".format(dst, sefd_loc[0].split("/")[-1]))
                """
print("All missing data obsids")
obs_with_missing_profiles = list(dict.fromkeys(obs_with_missing_profiles))
obs_with_missing_profiles.sort()
print(obs_with_missing_profiles)

"""
print("Obs with no dpp directory")
obs_no_dpp = []
obs_wth_dpp = []
for ob in obs_with_missing_profiles:
    if os.path.isdir("/astro/mwavcs/vcs/{}/dpp".format(ob)):
        obs_wth_dpp.append(ob)
    else:
        obs_no_dpp.append(ob)
print(obs_no_dpp)
print("Obs with dpp directory")
print(obs_wth_dpp)
"""

