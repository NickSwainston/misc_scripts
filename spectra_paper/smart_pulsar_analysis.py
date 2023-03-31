import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math

from pulsar_spectra.spectral_fit import find_best_spectral_fit
from pulsar_spectra.catalogue import collect_catalogue_fluxes

df = pd.read_csv("{}/SMART_pulsars_flux_update.csv".format(os.path.dirname(os.path.realpath(__file__))))
df = df.sort_values(by='Probability Best', ascending=False)
print(df)

# grab model specific data frames
spl_df  = df[df["Model"] == "simple_power_law"]
bpl_df  = df[df["Model"] == "broken_power_law"]
hfto_df = df[df["Model"] == "high_frequency_cut_off_power_law"]
lfto_df = df[df["Model"] == "low_frequency_turn_over_power_law"]
dtos_df = df[df["Model"] == "double_turn_over_spectrum"]
none_df = df[df["Model"].isnull()]
pre_spl_df  = df[df["Model before MWA"] == "simple_power_law"]
pre_bpl_df  = df[df["Model before MWA"] == "broken_power_law"]
pre_hfto_df = df[df["Model before MWA"] == "high_frequency_cut_off_power_law"]
pre_lfto_df = df[df["Model before MWA"] == "low_frequency_turn_over_power_law"]
pre_dtos_df = df[df["Model before MWA"] == "double_turn_over_spectrum"]
pre_none_df = df[df["Model before MWA"].isnull()]

# See if there any model dependent correlations
sns.set_theme(style="ticks")
f, ax = plt.subplots(figsize=(10, 5))
ax.xaxis.grid(True)
#ax.set_xscale("log")
# for column in ("N data flux", "ATNF Period (s)", "ATNF DM", "ATNF B_surf (G)"):
#     plt.clf()
#     sns.boxplot(x=column, y="Model", data=df,
#         order=[
#             "simple_power_law",
#             "broken_power_law",
#             "log_parabolic_spectrum",
#             "high_frequency_cut_off_power_law",
#             "low_frequency_turn_over_power_law",
#         ],
#         whis=[0, 100], width=.6, palette="vlag")

#     # Add in points to show each observation
#     sns.stripplot(x=column, y="Model", data=df,
#         order=[
#             "simple_power_law",
#             "broken_power_law",
#             "log_parabolic_spectrum",
#             "high_frequency_cut_off_power_law",
#             "low_frequency_turn_over_power_law",
#         ],
#         size=4, color=".3", linewidth=0)

#     # Tweak the visual presentation
#     #ax.set(ylabel="")
#     #plt.xticks(list(np.arange(0,80,10)))
#     sns.despine(trim=True, left=True)
#     plt.tight_layout(pad=0.5)
#     plt.savefig(f'model_{column.replace(" ", "_")}.png', bbox_inches='tight')


print("Before MWA")
print(f"{len(spl_df)} simple power law pulsars")
print(f"Average spectral index: {spl_df['pl_a'].mean()}")
print("\nWith MWA")
print(f"{len(pre_spl_df)} simple power law pulsars")
print(f"Average spectral index: {pre_spl_df['pre_pl_a'].mean()}")



# Count model distributions
print("\nPulsars for each model")
print("Before MWA\n----------")
print(f"PL:  {len(pre_spl_df )}")
print(f"BPL: {len(pre_bpl_df )}")
print(f"HFCO:{len(pre_hfto_df)}")
print(f"LFTO:{len(pre_lfto_df)}")
print(f"DTOS:{len(pre_dtos_df)}")
print(f"None:{len(pre_none_df)}")
print("\nWith MWA\n----------")
print(f"PL:  {len(spl_df )}")
print(f"BPL: {len(bpl_df )}")
print(f"HFCO:{len(hfto_df)}")
print(f"LFTO:{len(lfto_df)}")
print(f"DTOS:{len(dtos_df)}")
print(f"None:{len(none_df)}")
# print("\nBest fit pulsars for each model")
# print(f"PL:  {pulsars_pl[:4]}")
# print(f"     {pulsars_pl[4:8]}")
# print(f"BPL: {pulsars_bpl[:4]}")
# print(f"     {pulsars_bpl[4:8]}")
# print(f"LPS: {pulsars_lps[:4]}")
# print(f"     {pulsars_lps[4:8]}")
# print(f"HFCO:{pulsars_hfco[:4]}")
# print(f"     {pulsars_hfco[4:8]}")
# print(f"LFTO:{pulsars_lfto[:4]}")
# print(f"     {pulsars_lfto[4:8]}")
# plot it
fig, ax = plt.subplots(figsize=(10,5))
X = np.arange(5)
ax.bar(X + 0.00, (len(pre_spl_df ), len(pre_bpl_df ), len(pre_hfto_df), len(pre_lfto_df), len(pre_dtos_df)), color = 'b', width = 0.25)
ax.bar(X + 0.25, (len(spl_df ), len(bpl_df ), len(hfto_df), len(lfto_df), len(dtos_df)), color = 'g', width = 0.25)
plt.xticks(
    X,
    (
        "simple\npower law",
        "broken\npower law",
        "high frequency\ncut off",
        "low frequency\nturn over",
        "double turn\nover spectrum",
    ),
    #rotation=-30
)
ax.set_ylabel("Number of pulsars")
ax.legend(labels=['Before SMART', 'Including SMART'])
plt.savefig("model_comparison.png", bbox_inches='tight')

nparams = {
    "simple_power_law": 2,
    "broken_power_law": 4,
    "high_frequency_cut_off_power_law": 3,
    "low_frequency_turn_over_power_law": 4,
    "double_turn_over_spectrum": 5,
}

# See what happends with modesl without low frequency data
pulsars_changed_new_low = []
pulsars_changed_more = []
nparam_change = {}
chang_sum = 0
nchanged = 0
for index, row in df.iterrows():
    # if row["Min freq before MWA (MHz)"] > 154.:
        # if row["Model before MWA"] == row["Model"] or (row.isnull()["Model before MWA"] and row.isnull()["Model"]):
        #     print(f'same: {row["Pulsar"]:10}: {row["Model before MWA"]:35}   -->  {row["Model"]:35}')
        # else:
        # print(row["Min freq before MWA (MHz)"])
    if row["Model before MWA"] != row["Model"] and not row.isnull()["Model before MWA"]:
        ndiff = nparams[row["Model"]] - nparams[row["Model before MWA"]]
        print(f'DIFF: {ndiff} {row["Pulsar"]:10}: {row["Model before MWA"]:35}  -->  {row["Model"]:35}')
        if row["Min freq before MWA (MHz)"] > 154.:
            pulsars_changed_new_low.append((row["Pulsar"], row))
        else:
            pulsars_changed_more.append((row["Pulsar"], row))

        if ndiff in nparam_change.keys():
            nparam_change[ndiff][0] += 1
            nparam_change[ndiff][1].append(f'{row["Pulsar"]}: {row["Model before MWA"]:35}  -->  {row["Model"]:35}')
        else:
            nparam_change[ndiff] = [1, [f'{row["Pulsar"]}: {row["Model before MWA"]:35}  -->  {row["Model"]:35}']]
        chang_sum += ndiff
        nchanged += 1
print("")
for ndiff in sorted(nparam_change.keys()):
    print(f"nparam diff: {ndiff} for {nparam_change[ndiff][0]} pulsars")
    for mod_change in nparam_change[ndiff][1]:
        print(f"    {mod_change}")
print(f"Average nparam diff: {chang_sum/nchanged:.2f}\n")




cat_list_no_smart = collect_catalogue_fluxes(exclude=["Bhat_2022"])
cat_list_all = collect_catalogue_fluxes()

# # plot changed models with new lows
# nrows = math.ceil(len(pulsars_changed_new_low)/2)
# fig, ax = plt.subplots(nrows=nrows, ncols=2, figsize=(4*2,3*nrows))
# i = 0
# for pulsar, row in pulsars_changed_new_low:
#     # print(f"\n{pulsar} {row['Model before MWA']} --> {row['Model']}")
#     ri = i//2
#     ci = i%2
#     ax[ri][ci].set_title(pulsar)#, pad=30)
#     # Before smart
#     freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_no_smart[pulsar]
#     pre_models, pre_iminuit_results, pre_fit_infos, pre_p_best, pre_p_catagory = find_best_spectral_fit(
#         pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
#         plot_best=True, alternate_style=True, axis=ax[ri][ci], secondary_fit=True)
#     # with smart
#     freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_all[pulsar]
#     models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(
#         pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
#         plot_best=True, alternate_style=True, axis=ax[ri][ci],
#         ref_markers={"Bhat_2022": ("k",       "d", 7*0.7)})
#     i += 1
# if len(pulsars_changed_new_low)%2 == 1:
#     ax[-1][1].axes.get_xaxis().set_visible(False)
#     ax[-1][1].axes.get_yaxis().set_visible(False)
#     ax[-1][1].set_frame_on(False)
# fig.tight_layout()
# plt.savefig(f"smart_model_changed_new_low.png", bbox_inches='tight', dpi=300)

# # Plot changed models which aren't new lows
# nrows = math.ceil(len(pulsars_changed_more)/2)
# fig, ax = plt.subplots(nrows=nrows, ncols=2, figsize=(4*2,3*nrows))
# i = 0
# for pulsar, row in pulsars_changed_more:
#     # print(f"\n{pulsar} {row['Model before MWA']} --> {row['Model']}")
#     ri = i//2
#     ci = i%2
#     ax[ri][ci].set_title(pulsar)#, pad=30)
#     # Before smart
#     freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_no_smart[pulsar]
#     pre_models, pre_iminuit_results, pre_fit_infos, pre_p_best, pre_p_catagory = find_best_spectral_fit(
#         pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
#         plot_best=True, alternate_style=True, axis=ax[ri][ci], secondary_fit=True)
#     # with smart
#     freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_all[pulsar]
#     models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(
#         pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
#         plot_best=True, alternate_style=True, axis=ax[ri][ci],
#         ref_markers={"Bhat_2022": ("k",       "d", 7*0.7)})
#     i += 1
# if len(pulsars_changed_more)%2 == 1:
#     ax[-1][1].axes.get_xaxis().set_visible(False)
#     ax[-1][1].axes.get_yaxis().set_visible(False)
#     ax[-1][1].set_frame_on(False)
# fig.tight_layout()
# plt.savefig(f"smart_model_changed_more.png", bbox_inches='tight', dpi=300)

# Plot all changed models
# all_changes = pulsars_changed_more + pulsars_changed_new_low
# all_changes.sort()
# nplots = math.ceil(len(all_changes)/2)
# nrows = 4
# i = 0
# for pulsar, row in all_changes:
#     iplot = i//8
#     print(f"{i:2d}: {pulsar} {iplot}")
#     if i%8 == 0:
#         fig, ax = plt.subplots(nrows=nrows, ncols=2, figsize=(4*2,3*nrows))
#     # print(f"\n{pulsar} {row['Model before MWA']} --> {row['Model']}")
#     ri = (i-iplot*8)//2
#     ci = i%2
#     ax[ri][ci].set_title(pulsar)#, pad=30)
#     # Before smart
#     freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_no_smart[pulsar]
#     pre_models, pre_iminuit_results, pre_fit_infos, pre_p_best, pre_p_catagory = find_best_spectral_fit(
#         pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
#         plot_best=True, alternate_style=True, axis=ax[ri][ci], secondary_fit=True)
#     # with smart
#     freq_all, band_all, flux_all, flux_err_all, ref_all = cat_list_all[pulsar]
#     models, iminuit_results, fit_infos, p_best, p_catagory = find_best_spectral_fit(
#         pulsar, freq_all, band_all, flux_all, flux_err_all, ref_all,
#         plot_best=True, alternate_style=True, axis=ax[ri][ci],
#         ref_markers={"Bhat_2022": ("k",       "d", 7*0.7)})
#     if i%8 == 7:
#         fig.tight_layout()
#         plt.savefig(f"smart_model_changed_{iplot}.png", bbox_inches='tight', dpi=300)
#     i += 1
# if len(all_changes)%2 == 1:
#     ax[-1][1].axes.get_xaxis().set_visible(False)
#     ax[-1][1].axes.get_yaxis().set_visible(False)
#     ax[-1][1].set_frame_on(False)
# fig.tight_layout()
# plt.savefig(f"smart_model_changed_{iplot}.png", bbox_inches='tight', dpi=300)

# Make a summary latex table
short_names = {
    "simple_power_law":"SPL",
    "broken_power_law":"BPL",
    "high_frequency_cut_off_power_law":"HFCO",
    "low_frequency_turn_over_power_law":"LFTO",
    "double_turn_over_spectrum":"DTOS",
}
df = df.sort_values(by='Pulsar')
print("\cline{1-5}  \cline{7-11}")
print(f'{"PSR Jname"} & {"N flux"} & {"Original"} & {"New"} & {"Min Freq"} & & {"PSR Jname"} & {"N flux"} & {"Original"} & {"New"} & {"Min Freq"} \\\\')
print(f' & & {"Model"} & {"Model"} & {"(MHz)"} & &  & & {"Model"} & {"Model"} & {"(MHz)"}\\\\')
print("\cline{1-5}  \cline{7-11}")
low_flux_sum = 0
for pair1, pair2 in zip(df[:60].iterrows(), df[60:].iterrows()):
    row1 = pair1[1]
    pulsar1 = row1["Pulsar"]
    if not pd.isna(row1["Model before MWA"]):
        pre_model1 = short_names[row1["Model before MWA"]]
    else:
        pre_model1 = "-"
    if not pd.isna(row1["Model"]):
        model1 = short_names[row1["Model"]]
    else:
        model1 = "-"
    min_freq1 = int(row1["Min freq before MWA (MHz)"])
    if min_freq1 > 156:
        min_freq1 = f"{{\\bf {min_freq1}}}"
    ndata1 = row1["N data flux"] - 1
    if pre_model1 != model1:
        models1 = f"{{\\bf {pre_model1}}} & {{\\bf {model1}}}"
    else:
        models1 = f"{pre_model1} & {model1}"

    row2 = pair2[1]
    pulsar2 = row2["Pulsar"]
    if not pd.isna(row2["Model before MWA"]):
        pre_model2 = short_names[row2["Model before MWA"]]
    else:
        pre_model2 = "-"
    if not pd.isna(row2["Model"]):
        model2 = short_names[row2["Model"]]
    else:
        model2 = "-"
    min_freq2 = int(row2["Min freq before MWA (MHz)"])
    if min_freq2 > 156:
        low_flux_sum += 1
        min_freq2 = f"{{\\bf {min_freq2}}}"
    ndata2 = row2["N data flux"] -1
    if pre_model2 != model2:
        models2 = f"{{\\bf {pre_model2}}} & {{\\bf {model2}}}"
    else:
        models2 = f"{pre_model2} & {model2}"

    print(f'{pulsar1} & {ndata1} & {models1} & {min_freq1} & & {pulsar2} & {ndata2} & {models2} & {min_freq2}\\\\')
print("\cline{1-5}  \cline{7-11}")
print(f"Pulsars that SMART was lowest flux desnity measurement: {low_flux_sum}")

