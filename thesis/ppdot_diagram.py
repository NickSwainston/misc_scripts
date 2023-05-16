from psrqpy import QueryATNF
from pulsar_spectra.catalogue import collect_catalogue_fluxes

cat_list = collect_catalogue_fluxes()

# make basic plot
query = QueryATNF()
fig, ax = query.ppdot(showtypes=['BINARY'])

# make a list of jankwoski pulsars and my pulsars
jankowski = []
mine = []
for pulsar in cat_list.keys():
    if "Jankowski_2018" in cat_list[pulsar][4]:
        jankowski.append(pulsar)
    if len(cat_list[pulsar][0]) > 3:
        mine.append(pulsar)

# grab the data for jankwoski pulsars
periods = []
pdots = []
for i, pulsar in enumerate(query["PSRJ"]):
    if pulsar in jankowski:
        # record
        periods.append(query["P0"][i])
        pdots.append(query["P1"][i])

# Plot over the pulsars from jankowski
ax.loglog(
    periods,
    pdots,
    label="Jankowski et al. (2018)",
    marker='.',
    linestyle='none',
)
ax.legend(loc='upper left')
fig.savefig("ppdot_jankowski.png", dpi=300)


# grab the data for my pulsars
periods = []
pdots = []
for i, pulsar in enumerate(query["PSRJ"]):
    if pulsar in mine:
        # record
        periods.append(query["P0"][i])
        pdots.append(query["P1"][i])

# Plot over the pulsars from jankowski
ax.loglog(
    periods,
    pdots,
    label="This work",
    marker='.',
    linestyle='none',
)
ax.legend(loc='upper left')
fig.savefig("ppdot_mine.png", dpi=300)
