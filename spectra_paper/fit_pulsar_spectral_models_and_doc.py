import pandas as pd
import os
import matplotlib.pyplot as plt
import psrqpy
import numpy as np
import shutil
from PIL import Image
import glob
import torch.multiprocessing as mp
from functools import partial
from tqdm import tqdm
import yaml

from pulsar_spectra.spectral_fit import find_best_spectral_fit, estimate_flux_density
from pulsar_spectra.catalogue import collect_catalogue_fluxes, CAT_DIR
from pulsar_spectra.models import model_settings

df = pd.read_csv("{}/../survey_paper/SMART_pulsars.csv".format(os.path.dirname(os.path.realpath(__file__))))
pulsar_obsid_df = pd.read_csv("{}/pulsar_best_obs.csv".format(os.path.dirname(os.path.realpath(__file__))))

#print(df['Pulsar'].tolist())
#pulsars = list(dict.fromkeys(df['Pulsar'].tolist()))
with open(f"{CAT_DIR}/Bhat_2022.yaml", "r") as stream:
    cat_dict = yaml.safe_load(stream)
pulsars = cat_dict.keys()
print(len(pulsars))

cat_list = collect_catalogue_fluxes(exclude=["Bhat_2022"])
query = psrqpy.QueryATNF().pandas

results_record = []

#for output csv
output_df = pd.DataFrame(
    columns=[
        "Pulsar",
        "ObservationID",
        "ATNF Period (s)",
        "ATNF DM",
        "ATNF B_surf (G)",
        "ATNF E_dot (ergs/s)",
        "Offset (degrees)",
        "Flux Density (mJy)",
        "Flux Density Uncertainty (mJy)",
        "Flux Density Scintilation Uncertainty (mJy)",
        "Estimated Flux Density (mJy)",
        "Estimated Flux Density Uncertainty (mJy)",
        "Model",
        "Model before MWA",
        "Probability Best",
        "Min freq before MWA (MHz)",
        "Max freq before MWA (MHz)",
        "N data flux",
        "pl_a"      ,
        "pl_u_a"      ,
        "pl_c"      ,
        "pl_u_c"      ,
        "bpl_vb"    ,
        "bpl_u_vb"    ,
        "bpl_a1"    ,
        "bpl_u_a1"    ,
        "bpl_a2"    ,
        "bpl_u_a2"    ,
        "bpl_c"     ,
        "bpl_u_c"     ,
        "hfco_vc"   ,
        "hfco_u_vc"   ,
        "hfco_c"    ,
        "hfco_u_c"    ,
        "lfto_vpeak"   ,
        "lfto_u_vpeak" ,
        "lfto_a"    ,
        "lfto_u_a"    ,
        "lfto_c"    ,
        "lfto_u_c"    ,
        "lfto_beta" ,
        "lfto_u_beta" ,
        "dtos_vpeak"   ,
        "dtos_u_vpeak" ,
        "dtos_vc"   ,
        "dtos_u_vc" ,
        "dtos_a"    ,
        "dtos_u_a"    ,
        "dtos_c"    ,
        "dtos_u_c"    ,
        "dtos_beta" ,
        "dtos_u_beta" ,
        "pre_pl_a"      ,
        "pre_pl_u_a"      ,
        "pre_pl_c"      ,
        "pre_pl_u_c"      ,
        "pre_bpl_vb"    ,
        "pre_bpl_u_vb"    ,
        "pre_bpl_a1"    ,
        "pre_bpl_u_a1"    ,
        "pre_bpl_a2"    ,
        "pre_bpl_u_a2"    ,
        "pre_bpl_c"     ,
        "pre_bpl_u_c"     ,
        "pre_hfco_vc"   ,
        "pre_hfco_u_vc"   ,
        "pre_hfco_c"    ,
        "pre_hfco_u_c"    ,
        "pre_lfto_vpeak"   ,
        "pre_lfto_u_vpeak"   ,
        "pre_lfto_a"    ,
        "pre_lfto_u_a"    ,
        "pre_lfto_c"    ,
        "pre_lfto_u_c"    ,
        "pre_lfto_beta" ,
        "pre_lfto_u_beta" ,
        "pre_dtos_vpeak"   ,
        "pre_dtos_u_vpeak" ,
        "pre_dtos_vc"   ,
        "pre_dtos_u_vc" ,
        "pre_dtos_a"    ,
        "pre_dtos_u_a"    ,
        "pre_dtos_c"    ,
        "pre_dtos_u_c"    ,
        "pre_dtos_beta" ,
        "pre_dtos_u_beta" ,
    ]
)

model_dict = model_settings()
def fit_and_plot(pulsar):
#for pulsar in pulsars:
    scale_figure = 0.9
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5*scale_figure,4*scale_figure))
    # if pulsar != "J1136+1551":
    #      continue
    # if os.path.exists(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/{pulsar}.rst"):
    #     continue
    # Grab all rows with pulsar then grab their fluxes
    pulsar_df = df.loc[df['Pulsar'] == pulsar]
    mwa_freqs = []
    mwa_fluxs = []
    mwa_flux_errors = []
    pulsar_plots = []
    models = None
    pre_models = None
    estimate_flux = ""
    estimate_flux_err = ""

    for index, row in pulsar_df.iterrows():
        if row["Plot location"] == "" or isinstance(row["Plot location"], float):
            continue
        if not np.isnan(row['Flux Density (mJy)']) and row['Flux Density (mJy)'] != 0:
            mwa_freqs.append(154.24)
            mwa_fluxs.append(row['Flux Density (mJy)'])
            mwa_flux_errors.append(row['Flux Density Uncertainty (mJy)'])

        # Move files
        on_pulse_fit = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/on_pulse_plots/{row['ObservationID']}_{pulsar}_*_bins_gaussian_components.png")
        if len(on_pulse_fit) == 0:
            on_pulse_fit = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../survey_paper/{row['ObservationID']}_{pulsar}_*_bins_gaussian_components.png")
            if len(on_pulse_fit) == 0:
                on_pulse_fit = [""]
            else:
                #os.rename(on_pulse_fit[0], f"{os.path.dirname(os.path.realpath(__file__))}/../docs/on_pulse_plots/{on_pulse_fit[0]}")
                shutil.copyfile(on_pulse_fit[0], f"{os.path.dirname(os.path.realpath(__file__))}/../docs/on_pulse_plots/{os.path.basename(on_pulse_fit[0])}")
                on_pulse_fit = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/on_pulse_plots/{row['ObservationID']}_{pulsar}_*_bins_gaussian_components.png")

        if row["Plot location"].endswith("ps"):
            basename = row["Plot location"].split("/")[-1][:-2]
            png_basename = f"{basename}png"
            detection_plot = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{png_basename}")
            if len(detection_plot) == 0:
                # doesn't exist so make a png
                print(f"gs -sDEVICE=eps2write -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -r600 -sDEVICE=pngalpha -sOutputFile={os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{png_basename} {row['Plot location']}")
                os.system(f"gs -sDEVICE=eps2write -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -r600 -sDEVICE=pngalpha -sOutputFile={os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{png_basename} {row['Plot location']}")
                img = Image.open(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{png_basename}")
                # rotate by 90 degrees
                rot_img = img.transpose(Image.ROTATE_270)
                rot_img.save(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{png_basename}")
                detection_plot = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{png_basename}")
        else:
            basename = os.path.basename(row["Plot location"])
            detection_plot = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{basename}")
            if len(detection_plot) == 0:
                # cp
                print(row["Plot location"])
                shutil.copyfile(row["Plot location"], f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{os.path.basename(row['Plot location'])}")
                detection_plot = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{basename}")
        pulsar_plots.append((detection_plot[0], on_pulse_fit[0]))

    freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list[pulsar]
    freqs = freq_all + [154.24]
    fit_range = (np.log10(min(freqs)), np.log10(max(freqs)))

    pre_pl_a      = None
    pre_pl_u_a      = None
    pre_pl_c      = None
    pre_pl_u_c      = None
    pre_bpl_vb    = None
    pre_bpl_u_vb    = None
    pre_bpl_a1    = None
    pre_bpl_u_a1    = None
    pre_bpl_a2    = None
    pre_bpl_u_a2    = None
    pre_bpl_c     = None
    pre_bpl_u_c     = None
    pre_hfco_vc   = None
    pre_hfco_u_vc   = None
    pre_hfco_c    = None
    pre_hfco_u_c    = None
    pre_lfto_vpeak   = None
    pre_lfto_u_vpeak   = None
    pre_lfto_a    = None
    pre_lfto_u_a    = None
    pre_lfto_c    = None
    pre_lfto_u_c    = None
    pre_lfto_beta = None
    pre_lfto_u_beta = None
    pre_dtos_vpeak   = None
    pre_dtos_u_vpeak   = None
    pre_dtos_vc   = None
    pre_dtos_u_vc   = None
    pre_dtos_a    = None
    pre_dtos_u_a    = None
    pre_dtos_c    = None
    pre_dtos_u_c    = None
    pre_dtos_beta = None
    pre_dtos_u_beta = None
    pl_a      = None
    pl_u_a      = None
    pl_c      = None
    pl_u_c      = None
    bpl_vb    = None
    bpl_u_vb    = None
    bpl_a1    = None
    bpl_u_a1    = None
    bpl_a2    = None
    bpl_u_a2    = None
    bpl_c     = None
    bpl_u_c     = None
    hfco_vc   = None
    hfco_u_vc   = None
    hfco_c    = None
    hfco_u_c    = None
    lfto_vpeak   = None
    lfto_u_vpeak   = None
    lfto_a    = None
    lfto_u_a    = None
    lfto_c    = None
    lfto_u_c    = None
    lfto_beta = None
    lfto_u_beta = None
    dtos_vpeak   = None
    dtos_u_vpeak   = None
    dtos_vc   = None
    dtos_u_vc   = None
    dtos_a    = None
    dtos_u_a    = None
    dtos_c    = None
    dtos_u_c    = None
    dtos_beta = None
    dtos_u_beta = None
    if len(freq_all) > 0:
        pre_models, pre_iminuit_results, pre_fit_infos, pre_p_best, pre_p_catagory = find_best_spectral_fit(
            pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
            plot_best=True, alternate_style=True, axis=ax, secondary_fit=True
        )
    else:
        pre_models = pre_iminuit_results = pre_fit_infos = pre_p_best = pre_p_catagory = None
    if pre_models is not None:
        estimate_flux, estimate_flux_err = estimate_flux_density(154.24, pre_models, pre_iminuit_results)

        # record model specific bits
        if pre_models == "simple_power_law":
            pre_pl_a = pre_iminuit_results.values["a"]
            pre_pl_u_a = pre_iminuit_results.errors["a"]
            pre_pl_c = pre_iminuit_results.values["c"]
            pre_pl_u_c = pre_iminuit_results.errors["c"]
        elif pre_models == "broken_power_law":
            #vb, a1, a2, b
            pre_bpl_vb = pre_iminuit_results.values["vb"]
            pre_bpl_u_vb = pre_iminuit_results.errors["vb"]
            pre_bpl_a1 = pre_iminuit_results.values["a1"]
            pre_bpl_u_a1 = pre_iminuit_results.errors["a1"]
            pre_bpl_a2 = pre_iminuit_results.values["a2"]
            pre_bpl_u_a2 = pre_iminuit_results.errors["a2"]
            pre_bpl_c = pre_iminuit_results.values["c"]
            pre_bpl_u_c = pre_iminuit_results.errors["c"]
        elif pre_models == "high_frequency_cut_off_power_law":
            pre_hfco_vc = pre_iminuit_results.values["vc"]
            pre_hfco_u_vc = pre_iminuit_results.errors["vc"]
            pre_hfco_c = pre_iminuit_results.values["c"]
            pre_hfco_u_c = pre_iminuit_results.errors["c"]
        elif pre_models == "low_frequency_turn_over_power_law":
            #  vc, a, b, beta
            pre_lfto_vpeak = pre_iminuit_results.values["vpeak"]
            pre_lfto_u_vpeak = pre_iminuit_results.errors["vpeak"]
            pre_lfto_a = pre_iminuit_results.values["a"]
            pre_lfto_u_a = pre_iminuit_results.errors["a"]
            pre_lfto_c = pre_iminuit_results.values["c"]
            pre_lfto_u_c = pre_iminuit_results.errors["c"]
            pre_lfto_beta = pre_iminuit_results.values["beta"]
            pre_lfto_u_beta = pre_iminuit_results.errors["beta"]
        elif pre_models == "double_turn_over_spectrum":
            #  vc, a, b, beta
            pre_dtos_vc = pre_iminuit_results.values["vc"]
            pre_dtos_u_vc = pre_iminuit_results.errors["vc"]
            pre_dtos_vc = pre_iminuit_results.values["vpeak"]
            pre_dtos_u_vc = pre_iminuit_results.errors["vpeak"]
            pre_dtos_a = pre_iminuit_results.values["a"]
            pre_dtos_u_a = pre_iminuit_results.errors["a"]
            pre_dtos_c = pre_iminuit_results.values["c"]
            pre_dtos_u_c = pre_iminuit_results.errors["c"]
            pre_dtos_beta = pre_iminuit_results.values["beta"]
            pre_dtos_u_beta = pre_iminuit_results.errors["beta"]


    # calc offset
    # find obsid using for this pulsar
    #pulsar_obsid_df = pulsar_obsid_df[pulsar_obsid_df['Jname'].str.contains(pulsar)]
    this_df = pulsar_obsid_df.loc[pulsar == pulsar_obsid_df['Jname']].reset_index()
    obsid = this_df['ObsID'][0]
    query_id = list(query['PSRJ']).index(pulsar)
    # obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels = get_common_obs_metadata(obsid)
    # print("coords")
    # print(ra,dec)
    # obs_beam = SkyCoord(ra, dec, unit=(u.deg,u.deg))
    # query_id = list(query['PSRJ']).index(pulsar)
    # print(query["RAJ"][query_id], query["DECJ"][query_id])
    # pulsar_coord = SkyCoord(query["RAJ"][query_id], query["DECJ"][query_id], unit=(u.hourangle,u.deg))
    # offset = pulsar_coord.separation(obs_beam).deg
    offset = None


    print(f"\n{pulsar}")
    if len(mwa_fluxs) == 0:
        print(f"No fluxes")
    else:

        # Adjust uncertanty to take into account scintillation and number of detections
        # Average data
        S_mean = np.mean(mwa_fluxs)
        u_S_mean = np.sqrt(np.sum(np.array(mwa_flux_errors)**2))

        # using 728 MHz values from table 4
        a = -0.47
        b = 0.21
        d0 = 200
        d = float(row['ATNF DM'])
        # Equation 18 modultaion index
        m_r_v = b* (d/d0)**a

        # Equation 4
        u_scint = m_r_v * S_mean

        # Equation 2 robust standard deviation computed using the interquartile range
        std_r_v = 0.9183 * (np.quantile(mwa_fluxs, 0.75) - np.quantile(mwa_fluxs, 0.25))

        N = len(mwa_fluxs)
        # Equation 3
        u_S = np.sqrt( u_S_mean**2 + std_r_v**2/N + (6/(5*N) - 1/5)*u_scint**2)

        #print(u_S, mwa_flux_errors)

        # freq_all     = np.array(mwa_freqs                + cat_list[pulsar][0])
        # flux_all     = np.array(mwa_fluxs                + cat_list[pulsar][1])
        # flux_err_all = np.array(mwa_flux_errors          + cat_list[pulsar][2])
        # ref_all      = np.array(["SMART"]*len(mwa_freqs) + cat_list[pulsar][3])
        freq_all     = np.array([154.24]  + cat_list[pulsar][0])
        band_all     = np.array([30.72]   + cat_list[pulsar][1])
        flux_all     = np.array([S_mean]  + cat_list[pulsar][2])
        flux_err_all = np.array([u_S]     + cat_list[pulsar][3])
        ref_all      = np.array(["SMART"] + cat_list[pulsar][4])
        #for freq, flux, flux_err in zip(freq_all, flux_all, flux_err_all):
            #print(freq, flux, flux_err)
        models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(
            pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
            plot_best=True, alternate_style=True, axis=ax
        )
        plt.tight_layout(pad=2.5)
        #plt.savefig(f"{pulsar}_fit.pdf", bbox_inches='tight', dpi=300)
        plt.savefig(f"{pulsar}_fit.png", bbox_inches='tight', dpi=300)
        models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all, plot_compare=True)

        if models is not None:
            if len(models) > 0:
                shutil.move(f"{pulsar}_fit.png", f"{os.path.dirname(os.path.realpath(__file__))}/../docs/best_fits/{pulsar}_fit.png")
                print(f"{pulsar}_comparison_fit.png",  f"{os.path.dirname(os.path.realpath(__file__))}/../docs/comparison_fits/{pulsar}_comparison_fit.png")
                shutil.move(f"{pulsar}_comparison_fit.png",  f"{os.path.dirname(os.path.realpath(__file__))}/../docs/comparison_fits/{pulsar}_comparison_fit.png")

                # Record data
                results_record.append((pulsar, row['ATNF DM'], models, iminuit_results, fit_infos, p_best, p_catagory, len(mwa_fluxs), S_mean, u_S, u_S_mean, u_scint, m_r_v))


                # record model specific bits
                if models == "simple_power_law":
                    pl_a = iminuit_results.values["a"]
                    pl_u_a = iminuit_results.errors["a"]
                    pl_c = iminuit_results.values["c"]
                    pl_u_c = iminuit_results.errors["c"]
                elif models == "broken_power_law":
                    #vb, a1, a2, b
                    bpl_vb = iminuit_results.values["vb"]
                    bpl_u_vb = iminuit_results.errors["vb"]
                    bpl_a1 = iminuit_results.values["a1"]
                    bpl_u_a1 = iminuit_results.errors["a1"]
                    bpl_a2 = iminuit_results.values["a2"]
                    bpl_u_a2 = iminuit_results.errors["a2"]
                    bpl_c = iminuit_results.values["c"]
                    bpl_u_c = iminuit_results.errors["c"]
                elif models == "high_frequency_cut_off_power_law":
                    hfco_vc = iminuit_results.values["vc"]
                    hfco_u_vc = iminuit_results.errors["vc"]
                    hfco_c = iminuit_results.values["c"]
                    hfco_u_c = iminuit_results.errors["c"]
                elif models == "low_frequency_turn_over_power_law":
                    #  vc, a, b, beta
                    lfto_vpeak = iminuit_results.values["vpeak"]
                    lfto_u_vpeak = iminuit_results.errors["vpeak"]
                    lfto_a = iminuit_results.values["a"]
                    lfto_u_a = iminuit_results.errors["a"]
                    lfto_c = iminuit_results.values["c"]
                    lfto_u_c = iminuit_results.errors["c"]
                    lfto_beta = iminuit_results.values["beta"]
                    lfto_u_beta = iminuit_results.errors["beta"]
                elif models == "double_turn_over_spectrum":
                    #  vc, a, b, beta
                    dtos_vc = iminuit_results.values["vc"]
                    dtos_u_vc = iminuit_results.errors["vc"]
                    dtos_vc = iminuit_results.values["vpeak"]
                    dtos_u_vc = iminuit_results.errors["vpeak"]
                    dtos_a = iminuit_results.values["a"]
                    dtos_u_a = iminuit_results.errors["a"]
                    dtos_c = iminuit_results.values["c"]
                    dtos_u_c = iminuit_results.errors["c"]
                    dtos_beta = iminuit_results.values["beta"]
                    dtos_u_beta = iminuit_results.errors["beta"]
        if len(cat_list[pulsar][0]) == 0:
            min_freq = None
            max_freq = None
        else:
            min_freq = min(cat_list[pulsar][0])
            max_freq = max(cat_list[pulsar][0])



        with open(f'{os.path.dirname(os.path.realpath(__file__))}/../docs/{pulsar}.rst', 'w') as file:
            file.write(f".. _{pulsar}:\n{pulsar}\n")
            file.write("="*len(pulsar))

            # Fit with out data
            if models is not None:
                print(iminuit_results.values)
                file.write(f'''

Best Fit
--------
.. image:: best_fits/{pulsar}_fit.png
  :width: 800

.. csv-table:: {pulsar} fit results
''')
                header_str = '   :header: "model",'
                data_str = f'   "{models}",'
                for p, v, e in zip(iminuit_results.parameters, iminuit_results.values, iminuit_results.errors):
                    if p.startswith('v'):
                        header_str += f'"{p} (MHz)",'
                        data_str += f'"{int(v/1e6):d}±{int(e/1e6):d}",'
                    else:
                        header_str += f'"{p}",'
                        data_str += f'"{v:.2f}±{e:.2f}",'
                file.write(f'''{header_str[:-1]}

{data_str[:-1]}''')
            else:
                file.write(f'''

Best Fit
--------
Only {len(mwa_freqs)} MWA data and {len(cat_list[pulsar][0])} cat data available
''')

            # Fit without our data
            if pre_models is not None:
                file.write(f'''

Fit Before MWA
--------------

.. csv-table:: {pulsar} before fit results
''')
                header_str = '   :header: "model",'
                data_str = f'   "{pre_models}",'
                for p, v, e in zip(pre_iminuit_results.parameters, pre_iminuit_results.values, pre_iminuit_results.errors):
                    if p.startswith('v'):
                        header_str += f'"{p} (MHz)",'
                        data_str += f'"{int(v/1e6):d}±{int(e/1e6):d}",'
                    else:
                        header_str += f'"{p}",'
                        data_str += f'"{v:.2f}±{e:.2f}",'
                file.write(f'''{header_str[:-1]}

{data_str[:-1]}''')
            file.write(f'''


Flux Density Results
--------------------
.. csv-table:: {pulsar} flux density total results
   :header: "N obs", "Flux Density (mJy)", "u_S_mean", "u_scint", "m_r_v"

   "{len(mwa_fluxs)}",  "{S_mean:.1f}±{u_S:.1f}", "{u_S_mean:.1f}", "{u_scint:.1f}", "{m_r_v:.3f}"

.. csv-table:: {pulsar} flux density individual results
   :header: "ObsID", "Flux Density (mJy)"

''')
            for index, row in pulsar_df.iterrows():
                file.write(f'''    "{row['ObservationID']}", "{row['Flux Density (mJy)']:.1f}±{row['Flux Density Uncertainty (mJy)']:.1f}"\n''')

            # Comparison fit
            if models is not None:
                file.write(f'''
Comparison Fit
--------------
.. image:: comparison_fits/{pulsar}_comparison_fit.png
  :width: 800
''')

            # Detection plots
            file.write(f'''
Detection Plots
---------------
''')
            for detection_plot, on_pulse_fit in pulsar_plots:
                file.write(f'''
.. image:: detection_plots/{os.path.basename(detection_plot)}
  :width: 800

.. image:: on_pulse_plots/{os.path.basename(on_pulse_fit)}
  :width: 800''')

        # Record data for csv
        #output_df = output_df.append({
        return {
            "Pulsar":pulsar,
            "ObservationID":obsid,
            "ATNF Period (s)": row['ATNF Period (s)'],
            "ATNF DM": row['ATNF DM'],
            "ATNF B_surf (G)":query["BSURF"][query_id],
            "ATNF E_dot (ergs/s)":query["EDOT"][query_id],
            "Offset (degrees)":offset,
            "Flux Density (mJy)":S_mean,
            "Flux Density Uncertainty (mJy)":u_S_mean,
            "Flux Density Scintilation Uncertainty (mJy)":u_S,
            "Estimated Flux Density (mJy)":estimate_flux,
            "Estimated Flux Density Uncertainty (mJy)":estimate_flux_err,
            "Model":models,
            "Model before MWA":pre_models,
            "Probability Best":p_best,
            "Min freq before MWA (MHz)":min_freq,
            "Max freq before MWA (MHz)":max_freq,
            "N data flux": len(flux_all),
            "pl_a"      : pl_a     ,
            "pl_u_a"      : pl_u_a     ,
            "pl_c"      : pl_c     ,
            "pl_u_c"      : pl_u_c     ,
            "bpl_vb"    : bpl_vb   ,
            "bpl_u_vb"    : bpl_u_vb   ,
            "bpl_a1"    : bpl_a1   ,
            "bpl_u_a1"    : bpl_u_a1   ,
            "bpl_a2"    : bpl_a2   ,
            "bpl_u_a2"    : bpl_u_a2   ,
            "bpl_c"     : bpl_c    ,
            "bpl_u_c"     : bpl_u_c    ,
            "hfco_vc"   : hfco_vc  ,
            "hfco_u_vc"   : hfco_u_vc  ,
            "hfco_c"    : hfco_c   ,
            "hfco_u_c"    : hfco_u_c   ,
            "lfto_vpeak"   : lfto_vpeak  ,
            "lfto_u_vpeak" : lfto_u_vpeak,
            "lfto_a"    : lfto_a   ,
            "lfto_u_a"    : lfto_u_a   ,
            "lfto_c"    : lfto_c   ,
            "lfto_u_c"    : lfto_u_c   ,
            "lfto_beta" : lfto_beta,
            "lfto_u_beta" : lfto_u_beta,
            "dtos_vpeak"  : dtos_vpeak  ,
            "dtos_u_vpeak": dtos_u_vpeak,
            "dtos_vc"     : dtos_vc  ,
            "dtos_u_vc"   : dtos_u_vc,
            "dtos_a"      : dtos_a   ,
            "dtos_u_a"    : dtos_u_a   ,
            "dtos_c"      : dtos_c   ,
            "dtos_u_c"    : dtos_u_c   ,
            "dtos_beta"   : dtos_beta,
            "dtos_u_beta" : dtos_u_beta,
            "pre_pl_a"      : pre_pl_a     ,
            "pre_pl_u_a"      : pre_pl_u_a     ,
            "pre_pl_c"      : pre_pl_c     ,
            "pre_pl_u_c"      : pre_pl_u_c     ,
            "pre_bpl_vb"    : pre_bpl_vb   ,
            "pre_bpl_u_vb"    : pre_bpl_u_vb   ,
            "pre_bpl_a1"    : pre_bpl_a1   ,
            "pre_bpl_u_a1"    : pre_bpl_u_a1   ,
            "pre_bpl_a2"    : pre_bpl_a2   ,
            "pre_bpl_u_a2"    : pre_bpl_u_a2   ,
            "pre_bpl_c"     : pre_bpl_c    ,
            "pre_bpl_u_c"     : pre_bpl_u_c    ,
            "pre_hfco_vc"   : pre_hfco_vc  ,
            "pre_hfco_u_vc"   : pre_hfco_u_vc  ,
            "pre_hfco_c"    : pre_hfco_c   ,
            "pre_hfco_u_c"    : pre_hfco_u_c   ,
            "pre_lfto_vpeak"   : pre_lfto_vpeak  ,
            "pre_lfto_u_vpeak"   : pre_lfto_u_vpeak  ,
            "pre_lfto_a"    : pre_lfto_a   ,
            "pre_lfto_u_a"    : pre_lfto_u_a   ,
            "pre_lfto_c"    : pre_lfto_c   ,
            "pre_lfto_u_c"    : pre_lfto_u_c   ,
            "pre_lfto_beta" : pre_lfto_beta,
            "pre_lfto_u_beta" : pre_lfto_u_beta,
            "pre_dtos_vpeak"  : pre_dtos_vpeak  ,
            "pre_dtos_u_vpeak": pre_dtos_u_vpeak,
            "pre_dtos_vc"     : pre_dtos_vc  ,
            "pre_dtos_u_vc"   : pre_dtos_u_vc,
            "pre_dtos_a"      : pre_dtos_a   ,
            "pre_dtos_u_a"    : pre_dtos_u_a   ,
            "pre_dtos_c"      : pre_dtos_c   ,
            "pre_dtos_u_c"    : pre_dtos_u_c   ,
            "pre_dtos_beta"   : pre_dtos_beta,
            "pre_dtos_u_beta" : pre_dtos_u_beta,
        }
        #, ignore_index=True)

pbar = tqdm(pulsars)
# freeze params/function as object
fc_ = partial(fit_and_plot)
# set number of processes
p = mp.Pool(8)
# runs mp with params on pbar
#results = p.imap(fc_, pbar)
results = list(p.imap(fc_, pbar))
#print(results)
# close out and join processes
p.close()
p.join()

output_df = pd.DataFrame(results)
output_df.to_csv('SMART_pulsars_flux_update.csv', index=False)


# Record summary results on homepage
with open(f'{os.path.dirname(os.path.realpath(__file__))}/../docs/index.rst', 'w') as file:
    file.write(f'''.. pulsar_spectra documentation master file, created by
   sphinx-quickstart on Sat Feb 12 11:03:47 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pulsar_spectra's documentation!
==========================================


Flux Density Results
--------------------
.. csv-table::
   :header: "Pulsar", "DM", "N obs", "Flux Density (mJy)", "u_S_mean", "u_scint", "m_r_v"

''')
    # setp up other tables while I'm looping
    simple_power_law = []
    broken_power_law = []
    log_parabolic_spectrum = []
    high_frequency_cut_off_power_law = []
    low_frequency_turn_over_power_law = []
    for pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory, n_obs, S_mean, u_S, u_S_mean, u_scint, m_r_v in results_record:
        file.write(f'   ":ref:`{pulsar}<{pulsar}>`", "{dm}", "{n_obs}",  "{S_mean:.1f}±{u_S:.1f}", "{u_S_mean:.1f}", "{u_scint:.1f}", "{m_r_v:.3f}"\n')

        #sort
        if models == "simple_power_law":
            simple_power_law.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models == "broken_power_law":
            broken_power_law.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models == "log_parabolic_spectrum":
            log_parabolic_spectrum.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models == "high_frequency_cut_off_power_law":
            high_frequency_cut_off_power_law.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models == "low_frequency_turn_over_power_law":
            low_frequency_turn_over_power_law.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))

    file.write(f'''
Single Power Law Results
------------------------
.. csv-table::
   :header: "Pulsar", "DM", "a", "b"

''')
    for pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory in simple_power_law:
        data_str = f'   ":ref:`{pulsar}<{pulsar}>`", "{dm}", '
        for p, v, e in zip(iminuit_results.parameters, iminuit_results.values, iminuit_results.errors):
            if p.startswith('v'):
                data_str += f'"{int(v/1e6):d}±{int(e/1e6):d}", '
            else:
                data_str += f'"{v:.2f}±{e:.2f}", '
        file.write(f'{data_str[:-2]}\n')

    file.write(f'''
Broken Power Law Results
------------------------
.. csv-table::
   :header: "Pulsar", "DM", "vb (MHz)", "a1", "a2", "b"

''')
    for pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory in broken_power_law:
        data_str = f'   ":ref:`{pulsar}<{pulsar}>`", "{dm}", '
        for p, v, e in zip(iminuit_results.parameters, iminuit_results.values, iminuit_results.errors):
            if p.startswith('v'):
                data_str += f'"{int(v/1e6):d}±{int(e/1e6):d}", '
            else:
                data_str += f'"{v:.2f}±{e:.2f}", '
        file.write(f'{data_str[:-2]}\n')

    file.write(f'''
Log-parabolic spectrum Results
------------------------------
.. csv-table::
   :header: "Pulsar", "DM", "a", "b", "c"

''')
    for pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory in log_parabolic_spectrum:
        data_str = f'   ":ref:`{pulsar}<{pulsar}>`", "{dm}", '
        for p, v, e in zip(iminuit_results.parameters, iminuit_results.values, iminuit_results.errors):
            if p.startswith('v'):
                data_str += f'"{int(v/1e6):d}±{int(e/1e6):d}", '
            else:
                data_str += f'"{v:.2f}±{e:.2f}", '
        file.write(f'{data_str[:-2]}\n')

    file.write(f'''
Power law with high-frequency cut-off Results
---------------------------------------------
.. csv-table::
   :header: "Pulsar", "DM", "vc (MHz)", "a", "b"

''')
    for pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory in high_frequency_cut_off_power_law:
        data_str = f'   ":ref:`{pulsar}<{pulsar}>`", "{dm}", '
        for p, v, e in zip(iminuit_results.parameters, iminuit_results.values, iminuit_results.errors):
            if p.startswith('v'):
                data_str += f'"{int(v/1e6):d}±{int(e/1e6):d}", '
            else:
                data_str += f'"{v:.2f}±{e:.2f}", '
        file.write(f'{data_str[:-2]}\n')

    file.write(f'''
Power law with low-frequency turn-over Results
----------------------------------------------
.. csv-table::
   :header: "Pulsar", "DM", "vc (MHz)", "a", "b", "beta"

''')
    for pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory in low_frequency_turn_over_power_law:
        data_str = f'   ":ref:`{pulsar}<{pulsar}>`", "{dm}", '
        for p, v, e in zip(iminuit_results.parameters, iminuit_results.values, iminuit_results.errors):
            if p.startswith('v'):
                data_str += f'"{int(v/1e6):d}±{int(e/1e6):d}", '
            else:
                data_str += f'"{v:.2f}±{e:.2f}", '
        file.write(f'{data_str[:-2]}\n')


    file.write(f'''
.. toctree::
   :maxdepth: 1
   :caption: Pulsar Fit Results:
   :glob:

   J*

''')
