import pandas as pd
import os
import matplotlib.pyplot as plt

df = pd.read_csv("{}/psr2_J0024-1932_follow_up.csv".format(os.path.dirname(os.path.realpath(__file__))))
obsids = []
sns = []
print(df.keys())
for i in range(38):
    if df["best pointing "][i] != "" and not pd.isnull(df["best pointing "][i]):
        print(df["best pointing "][i])
        print(df["SN"][i])
        obsids.append(int(df["obsid"][i]))
        sns.append(float(df["SN"][i]))

plt.scatter(obsids, sns)
#ax2.set_xticklabels((ax2.get_xticks()-59000).astype(int))
plt.xlabel('Time (GPS)')
plt.ylabel('S/N')
plt.savefig("psr_sn.png")
