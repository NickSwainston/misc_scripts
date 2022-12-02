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

cat_list_no_smart = collect_catalogue_fluxes(exclude=["Bhat_2022"])
cat_list_all = collect_catalogue_fluxes()
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
    print(f"\n{pulsar}")
    scale_figure = 0.9
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5*scale_figure,4*scale_figure))
    # if pulsar != "J1136+1551":
    #      continue
    # if os.path.exists(f"{os.path.dirname(os.path.realpath(__file__))}/../docs/{pulsar}.rst"):
    #     continue
    # Grab all rows with pulsar then grab their fluxes
    pulsar_df = df.loc[df['Pulsar'] == pulsar]

    # Set defaults
    models = None
    pre_models = None
    estimate_flux = ""
    estimate_flux_err = ""
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

    # Fit without SMART
    freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_no_smart[pulsar]
    pre_models, pre_iminuit_results, pre_fit_infos, pre_p_best, pre_p_catagory = find_best_spectral_fit(
        pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
        plot_best=True, alternate_style=True, axis=ax, secondary_fit=True
    )
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

    # Fit again with SMART
    freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_all[pulsar]
    models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(
        pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
        plot_best=True, alternate_style=True, axis=ax
    )
    plt.tight_layout(pad=2.5)
    #plt.savefig(f"{pulsar}_fit.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{pulsar}_fit.png", bbox_inches='tight', dpi=300)
    models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all, plot_compare=True)

    if models is not None:
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
    min_freq = min(cat_list_no_smart[pulsar][0])
    max_freq = max(cat_list_no_smart[pulsar][0])

    # Record data for csv
    #output_df = output_df.append({
    return {
        "Pulsar":pulsar,
        # "ObservationID":obsid,
        # "ATNF Period (s)": row['ATNF Period (s)'],
        # "ATNF DM": row['ATNF DM'],
        # "ATNF B_surf (G)":query["BSURF"][query_id],
        # "ATNF E_dot (ergs/s)":query["EDOT"][query_id],
        # "Offset (degrees)":offset,
        # "Flux Density (mJy)":S_mean,
        # "Flux Density Uncertainty (mJy)":u_S_mean,
        # "Flux Density Scintilation Uncertainty (mJy)":u_S,
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
