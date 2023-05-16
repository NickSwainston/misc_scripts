import psrqpy
import pandas
import numpy as np
import math
from scipy.stats import lognorm, norm
import seaborn as sns
import matplotlib.pyplot as plt

def fit_plot_distribution(dist, label):
    #print(dist[:10])
    norm_dist = np.log10(dist, dtype=np.float64)
    mu = np.mean(norm_dist)
    mu_exp = np.exp(mu)
    #print(mu, mu_exp, np.max(dist))
    #s1400 test
    t = np.linspace(np.min(dist), np.max(dist), 100)
    f, (ax1, ax2) = plt.subplots(1, 2, sharex='col', figsize=(10, 5))
    sns.histplot(dist, ax=ax1, kde=False)
    fitting_params_lognormal = lognorm.fit(dist, floc=0, scale=mu_exp)
    lognorm_dist_fitted = lognorm(*fitting_params_lognormal)
    print('Fitted Model X~LogNorm(mu={0:.1f}, sigma={1:.1f})'.format(lognorm_dist_fitted.mean(), lognorm_dist_fitted.std()))
    ax1.plot(t, lognorm_dist_fitted.pdf(t))

    #print(norm_dist[:10],norm_dist[-10:])
    t = np.linspace(np.min(norm_dist), np.max(norm_dist), 100)
    sns.histplot(norm_dist, ax=ax2, kde=False, stat='density')
    fitting_params_normal = norm.fit(norm_dist)
    norm_dist_fitted = norm(*fitting_params_normal)
    print('Fitted Model X~norm(mu={0:.1f}, sigma={1:.1f})'.format(norm_dist_fitted.mean(), norm_dist_fitted.std()))
    ax2.plot(t, norm_dist_fitted.pdf(t))

    #plt.show()
    plt.savefig("{}_dist_fit.png".format(label))

survey_list = ['MWA_SMART_shallow', 'MWA_SMART2020', 'LOFAR', 'PMSURV']
for survey in survey_list:
    x = np.genfromtxt("/home/nick/code/PsrPopPy/work_dir/{}_asci.txt".format(survey))
    print("\nFitting {} with {} pulsars".format(survey, len(x)))
    dist = x[1:,6]
    fit_plot_distribution(dist, survey)




print("\nRunning query")
query = psrqpy.QueryATNF().pandas
print("Finished query\n")

"""
headers = query.keys()
for h in headers:
    #if h[0] == 'S' and h[1].isdigit() and h[-1].isdigit():
    if h[0] == 'R':
        print(h)
"""


flux_queries = ["S40", "S50", "S60", "S80", 
                "S100", "S150", "S200", "S300", "S400",
                "S600", "S700", "S800", "S900",
                "S1400", "S1600", "S2000", "S3000",
                "S4000", "S5000", "S6000", "S8000"]

# Make a dictionary of each frequency with a list to be filled with luminosities
luminosity_results = {}
for f in flux_queries:
    luminosity_results[f] = []

# loop over pulsar and calc luminosity at each frequency
for pi in range(len(query)):
    dist = query['DIST'][pi] #kpc
    for f in flux_queries:
        flux = query[f][pi] #mJy
        """
        lum = query['R_LUM14'][pi]
        if not np.isnan(lum):
            if lum != 0.:
                luminosity_results[f].append(lum)
        """
        if not np.isnan(flux) and not np.isnan(dist):
            #print(flux, dist)
            lum = flux * math.pow(dist, 2)
            if not np.isfinite(lum):
                print(lum)
            if lum != 0.:
                luminosity_results[f].append(lum)

for f in flux_queries:
    if len(luminosity_results[f]) < 100:
        print("skipping {}: {}".format(f, len(luminosity_results[f])))
        continue
    else:
        print("\nFitting {} with {} pulsars".format(f, len(luminosity_results[f])))
    dist = luminosity_results[f]
    fit_plot_distribution(dist, f)
    