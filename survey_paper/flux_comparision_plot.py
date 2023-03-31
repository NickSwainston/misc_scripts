import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np


df = pd.read_csv("{}/SMART_pulsars_flux_update.csv".format(os.path.dirname(os.path.realpath(__file__))))
print(df)

measured_flux    = np.array(df['Flux Density (mJy)'].tolist())
measured_flux_u  = np.array(df['Flux Density Uncertainty (mJy)'].tolist())
predicted_flux   = np.array(df['Estimated Flux Density (mJy)'].tolist())
predicted_flux_u = np.array(df['Estimated Flux Density Uncertainty (mJy)'].tolist())
dms = df['ATNF DM'].tolist()
#antf_DM = df['ATNF DM'].tolist()

# Remove ones with any nans
mmf = []
mmfu = []
mpf = []
mpfu = []
for mf, mfu, pf, pfu in zip(measured_flux, measured_flux_u, predicted_flux, predicted_flux_u):
    if not (np.isnan(mf) or np.isnan(mfu) or np.isnan(pf) or np.isnan(pfu)):
        mmf.append(mf)
        mmfu.append(mfu)
        mpf.append(pf)
        mpfu.append(pfu)
masked_measured_flux    = np.array(mmf)
masked_measured_flux_u  = np.array(mmfu)
masked_predicted_flux   = np.array(mpf)
masked_predicted_flux_u = np.array(mpfu)
# Perfom comparison math
ratio = abs( masked_measured_flux / masked_predicted_flux )
print(ratio)
ratio_u = ratio * np.sqrt( (masked_predicted_flux_u/masked_predicted_flux)**2 + (masked_measured_flux_u/masked_measured_flux_u)**2 )

mean_ratio = np.mean(ratio)
mean_ratio_u = np.sqrt(np.sum(ratio_u**2))
print(f"Mean ratio = {mean_ratio:.2f} +/- {mean_ratio_u:.2f}")

# colour scattered pulsars
# colour = []
# for index, row in df.iterrows():
#     if row['Pulsar'] in ['J0820-4114', 'J0835-4510', 'J0907-5157', 'J1121-5444', 'J1202-5820']:
#         colour.append('r')
#     else:
#         colour.append('b')

print(min(measured_flux))
print(min(predicted_flux))
print(max(measured_flux))
print(max(predicted_flux))
print(len(predicted_flux), len(measured_flux), len(predicted_flux_u), len(measured_flux_u))

fig, ax = plt.subplots(1,1)

# plot expected line
max_sn = max(measured_flux + predicted_flux)
#min_sn = min(measured_flux + predicted_flux)
min_sn = 10
print([min_sn,max_sn], [min_sn, max_sn])

# Plot error bars
markersize = 2
makerwidth = 1
capsize = 3
fig, ax = plt.subplots()
ax.plot([min_sn, max_sn], [min_sn, max_sn], color='r')
erf = 2.
ax.plot([min_sn*erf, max_sn*erf], [min_sn, max_sn], color='orange')
ax.plot([min_sn, max_sn], [min_sn*erf, max_sn*erf], color='orange')
#else:
plt.scatter(
    predicted_flux,
    measured_flux,
    s=5,
    c=dms,
    cmap='cool',
)
ax.plot([min_sn, max_sn], [min_sn, max_sn], color='r')

plt.xlabel(r"Expected Flux Density, $S_{exp}$ (mJy)")
plt.ylabel(r"Measured Flux Density, $S_{meas}$ (mJy)")
plt.xlim([min_sn,max_sn])
plt.ylim([min_sn,max_sn])
plt.xscale("log")
plt.yscale("log")
ax.get_xaxis().set_major_formatter(ScalarFormatter())
ax.get_yaxis().set_major_formatter(ScalarFormatter())
plt.colorbar(label="DM")
plt.savefig("flux_comparison.png", bbox_inches='tight', dpi=1000)


for pf, mf, pfu, mfu, dm in zip(predicted_flux, measured_flux, predicted_flux_u, measured_flux_u, dms):
    (_, caps, _) = ax.errorbar(pf, mf,
                                xerr=pfu, yerr=mfu,
                                #c=antf_DM, cmap='cool',
                                color='black',
                                alpha=0.7,
                                fmt="o", markersize=markersize, capsize=capsize, zorder=0)
    for cap in caps:
        cap.set_markeredgewidth(makerwidth)
plt.scatter(
    predicted_flux,
    measured_flux,
    s=5,
    c=dms,
    cmap='cool',
)
plt.savefig("flux_comparison_error_bars.png", bbox_inches='tight', dpi=1000)

