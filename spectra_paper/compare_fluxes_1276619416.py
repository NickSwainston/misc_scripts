#! /usr/bin/env python

import os
import numpy as np
import csv

from pulsar_spectra import catalogues
from pulsar_spectra.spectral_fit import find_best_spectral_fit
from pulsar_spectra.catalogues import collect_catalogue_fluxes

import logging
logger = logging.getLogger(__name__)


cat_dict, cat_list = collect_catalogue_fluxes()
# 1276619416 results
#               jname, imaging flux (mJy), my flux
pulsars_data = [('J1820-0427', 545.428),
                ('J1825-0935', 176.222),
                ('J1834-0010', 107.99),
                ('J1834-0426', 259.982),
                ('J1849-0636', 84.704)]
for pulsar, image_flux in pulsars_data:
    print(f"\nFitting {pulsar}")
    
    # add data to dict but not list so it doesn't take it into account for fitting
    cat_dict[pulsar]["MWA imaging"] = {"Frequency MHz":[185],
                                     "Flux Density mJy":[image_flux],
                                     "Flux Density error mJy":[image_flux*.1]}

    freq_all = np.array(cat_list[pulsar][0])*1e6
    flux_all = np.array(cat_list[pulsar][1])*1e-3
    flux_err_all = np.array(cat_list[pulsar][2])*1e-3
    #print(freq_all, flux_all, flux_err_all)
    models, fit_results = find_best_spectral_fit(pulsar, freq_all, flux_all, flux_err_all, plot_best=True, data_dict=cat_dict[pulsar])
    print(models)