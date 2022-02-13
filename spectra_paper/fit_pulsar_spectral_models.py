import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import psrqpy
import numpy as np
import shutil
from PIL import Image
import glob

from pulsar_spectra.spectral_fit import find_best_spectral_fit
from pulsar_spectra.catalogues import collect_catalogue_fluxes

df = pd.read_csv("{}/../survey_paper/SMART_pulsars.csv".format(os.path.dirname(os.path.realpath(__file__))))

print(df['Pulsar'].tolist())
pulsars = list(dict.fromkeys(df['Pulsar'].tolist()))

cat_list = collect_catalogue_fluxes(exclude=["Xue_2017", "Bondonneau_2020", "Johnston_2021"])

results_record = []

for pulsar in pulsars:
    #if pulsar != "J0528+2200":
    # if pulsar != "J0030+0451":
    #    continue
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
            print(detection_plot)
            if len(detection_plot) == 0:
                # cp
                print(row["Plot location"])
                shutil.copyfile(row["Plot location"], f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{os.path.basename(row['Plot location'])}")
                detection_plot = glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/detection_plots/{basename}")
        pulsar_plots.append((detection_plot[0], on_pulse_fit[0]))

    freq_all     = np.array(cat_list[pulsar][0])
    flux_all     = np.array(cat_list[pulsar][1])
    flux_err_all = np.array(cat_list[pulsar][2])
    ref_all      = np.array(cat_list[pulsar][3])
    if len(freq_all) > 0:
        pre_models, pre_iminuit_results, pre_fit_infos, pre_p_best, pre_p_catagory = find_best_spectral_fit(pulsar, freq_all, flux_all, flux_err_all, ref_all, plot_compare=True, plot_best=True, alternate_style=True)
    else:
        pre_models = pre_iminuit_results = pre_fit_infos = pre_p_best = pre_p_catagory = None
    if pre_models is not None:
        if len(pre_models) > 0:
            shutil.move(f"{pulsar}_{pre_models[1]}_fit.png", f"{os.path.dirname(os.path.realpath(__file__))}/../docs/before_mwa/{pulsar}_{pre_models[1]}_fit.png")

    print(f"\n{pulsar}")
    if len(mwa_fluxs) == 0:
        print(f"No fluxes")
    else:

        # Adjust uncertanty to take into account scintillation and number of detections
        # Average data
        S_mean = np.mean(mwa_fluxs)
        u_S_mean = np.sqrt(np.sum(np.array(mwa_flux_errors)**2))
        # robust standard deviation com- puted using the interquartile range
        std_r_v = 0.9183 * (np.quantile(mwa_fluxs, 0.75) - np.quantile(mwa_fluxs, 0.25))
        # modultaion index
        #m_r = std_r_v / np.median(mwa_fluxs)
        # using 728 MHz values from table 4
        a = -0.47
        b = 0.21
        d0 = 200
        d = float(row['ATNF DM'])
        # Equation 18
        m_r_v = b* (d/d0)**a
        u_scint = m_r_v * S_mean

        N = len(mwa_fluxs)
        u_S = np.sqrt( u_S_mean**2 + std_r_v**2/N + (6/(5*N) - 1/5)*u_scint**2)
        #print(u_S, mwa_flux_errors)

        # freq_all     = np.array(mwa_freqs                + cat_list[pulsar][0])
        # flux_all     = np.array(mwa_fluxs                + cat_list[pulsar][1])
        # flux_err_all = np.array(mwa_flux_errors          + cat_list[pulsar][2])
        # ref_all      = np.array(["SMART"]*len(mwa_freqs) + cat_list[pulsar][3])
        freq_all     = np.array([154.24]  + cat_list[pulsar][0])
        flux_all     = np.array([S_mean]  + cat_list[pulsar][1])
        flux_err_all = np.array([u_S]     + cat_list[pulsar][2])
        ref_all      = np.array(["SMART"] + cat_list[pulsar][3])
        models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(pulsar, freq_all, flux_all, flux_err_all, ref_all, plot_compare=True, plot_best=True, alternate_style=True)
        if models is not None:
            if len(models) > 0:
                shutil.move(f"{pulsar}_{models[1]}_fit.png", f"{os.path.dirname(os.path.realpath(__file__))}/../docs/best_fits/{pulsar}_{models[1]}_fit.png")
                shutil.move(f"{pulsar}_comparison_fit.png",  f"{os.path.dirname(os.path.realpath(__file__))}/../docs/comparison_fits/{pulsar}_comparison_fit.png")

                # Record data
                results_record.append((pulsar, row['ATNF DM'], models, iminuit_results, fit_infos, p_best, p_catagory, len(mwa_fluxs), S_mean, u_S, u_S_mean, u_scint, m_r_v))


        with open(f'{os.path.dirname(os.path.realpath(__file__))}/../docs/{pulsar}.rst', 'w') as file:
            file.write(f".. _{pulsar}:\n{pulsar}\n")
            file.write("="*len(pulsar))

            # Fit with out data
            if models is not None:
                print(iminuit_results.values)
                file.write(f'''

Best Fit
--------
.. image:: best_fits/{pulsar}_{models[1]}_fit.png
  :width: 800

.. csv-table:: {pulsar} fit results
''')
                header_str = '   :header: "model",'
                data_str = f'   "{models[1]}",'
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
.. image:: before_mwa/{pulsar}_{pre_models[1]}_fit.png
  :width: 800

.. csv-table:: {pulsar} before fit results
''')
                header_str = '   :header: "model",'
                data_str = f'   "{pre_models[1]}",'
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
        if models[1] == "simple_power_law":
            simple_power_law.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models[1] == "broken_power_law":
            broken_power_law.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models[1] == "log_parabolic_spectrum":
            log_parabolic_spectrum.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models[1] == "high_frequency_cut_off_power_law":
            high_frequency_cut_off_power_law.append((pulsar, dm, models, iminuit_results, fit_infos, p_best, p_catagory))
        if models[1] == "low_frequency_turn_over_power_law":
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