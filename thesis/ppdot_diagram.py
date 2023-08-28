from psrqpy import QueryATNF

# make basic plot
query = QueryATNF()
fig, ax = query.ppdot(showtypes=['BINARY'])
beam = [
    "J0514-4408",
    "J0600-5756",
    "J0636-4549",
    "J0729-1448",
    "J0729-1836",
    "J0842-4851",
    "J0902-6325",
    "J1012-2337",
    "J1018-1642",
    "J1034-3224",
    "J1123-4844",
    "J1146-6030",
    "J1225-5556",
    "J1240-4124",
    "J1312-5402",
    "J1320-5359",
    "J2324-6054",
]

# grab the data for my pulsars
periods = []
pdots = []
for i, pulsar in enumerate(query["PSRJ"]):
    if pulsar in beam:
        # record
        periods.append(query["P0"][i])
        pdots.append(query["P1"][i])

# Plot over the pulsars from
ax.loglog(
    periods,
    pdots,
    label="Beamforming",
    marker='.',
    linestyle='none',
    markersize=15,
    markerfacecolor='orange',
)


image = [
    "J1231-6303",
    "J1527-5552",
    "J1600-5044",
    "J1614-5048",
    "J1623-4256",
    "J1633-5015",
    "J1639-4604",
    "J1644-4559",
    "J1721-3532",
    "J1807-2715",
    "J1823-1115",
    "J1827-0958",
    "J1832-0827",
    "J1833-0559",
    "J1836-1008",
    "J1903+0135",
    "J1751-2737",
    "J1748-3009",
]

# grab the data for my pulsars
periods = []
pdots = []
for i, pulsar in enumerate(query["PSRJ"]):
    if pulsar in image:
        # record
        periods.append(query["P0"][i])
        pdots.append(query["P1"][i])

# Plot over the pulsars from
ax.loglog(
    periods,
    pdots,
    label="Imaging",
    marker='.',
    linestyle='none',
    markersize=15,
    markerfacecolor='yellow',
)


both = [
    "J0737-3039A",
    "J0742-2822",
    "J0820-1350",
    "J0820-4114",
    "J0837+0610",
    "J0837-4135",
    "J0856-6137",
    "J0907-5157",
    "J0908-1739",
    "J0922+0638",
    "J0924-5302",
    "J1041-1942",
    "J1057-5226",
    "J1059-5742",
    "J1116-4122",
    "J1121-5444",
    "J1136+1551",
    "J1136-5525",
    "J1141-6545",
    "J1202-5820",
    "J1224-6407",
    "J1312-5402",
    "J1430-6623",
    "J1453-6413",
    "J1456-6843",
    "J1507-4352",
    "J1534-5334",
    "J1544-5308",
    "J1607-0032",
    "J1645-0317",
    "J1709-1640",
    "J1731-4744",
    "J1751-4657",
    "J1752-2806",
    "J1820-0427",
    "J1823-3106",
    "J1825-0935",
    "J1900-2600",
    "J1902-5105",
    "J1913-0440",
    "J1921+2153",
    "J1932+1059",
    "J1935+1616",
    "J1943-1237",
    "J2018+2839",
    "J2022+2854",
    "J2046-0421",
    "J2048-1616",
    "J2145-0750",
]

# grab the data for my pulsars
periods = []
pdots = []
for i, pulsar in enumerate(query["PSRJ"]):
    if pulsar in both:
        # record
        periods.append(query["P0"][i])
        pdots.append(query["P1"][i])

# Plot over the pulsars from
ax.loglog(
    periods,
    pdots,
    label="Both",
    marker='.',
    linestyle='none',
    markersize=15,
    markerfacecolor='purple',
)
ax.legend(loc='upper left')
fig.savefig("ppdot_mine.png", dpi=300)
