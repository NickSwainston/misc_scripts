import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


df = pd.read_csv("{}/SMART_pulsars.csv".format(os.path.dirname(os.path.realpath(__file__))))
print(df)

measured_sn = df['SN'].tolist()
measured_sn_u = df['SN Uncertainty'].tolist()
predicted_sn = df['Predicted SN'].tolist()
predicted_sn_u = df['Predicted SN Uncertainty'].tolist()
antf_DM = df['ATNF DM'].tolist()

print(min(measured_sn))
print(min(predicted_sn))

fig, ax = plt.subplots(1,1)

# plot expected line
max_sn = max(measured_sn + predicted_sn)
min_sn = min(measured_sn + predicted_sn)
plt.plot([min_sn,max_sn], [min_sn, max_sn], color='r')

# Plot error bars
markersize = 3
makerwidth = 1
capsize = 3

error_bars = False
if error_bars:
    (_, caps, _) = plt.errorbar(predicted_sn, measured_sn,
                                xerr=predicted_sn_u, yerr=measured_sn_u,
                                #c=antf_DM, cmap='cool',
                                fmt="o", markersize=markersize, capsize=capsize, zorder=0)
    for cap in caps:
        cap.set_markeredgewidth(makerwidth)
#else:
plt.scatter(predicted_sn, measured_sn, s=5,
            c=antf_DM, cmap='cool'
            )
plt.colorbar(spacing='uniform',
             label=r"ANTF DM (pc cm$^{-3}$)")

plt.xlabel("Predicted S/N")
plt.ylabel("Measured S/N")
plt.xlim([min_sn,max_sn])
plt.ylim([min_sn,max_sn])
plt.xscale("log")
plt.yscale("log")
ax.get_xaxis().set_major_formatter(ScalarFormatter())
ax.get_yaxis().set_major_formatter(ScalarFormatter())
plt.savefig("sn_comparison.png", bbox_inches='tight', dpi=1000)

