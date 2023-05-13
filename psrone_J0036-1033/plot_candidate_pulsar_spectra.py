import numpy as np
from pulsar_spectra.spectral_fit import find_best_spectral_fit
import pandas as pd

def average_flux(flux, flux_u):
    flux   = np.array(flux)
    flux_u = np.array(flux_u)

    fsum  = np.sum(flux)
    fusum = np.sqrt(np.sum(flux_u**2))

    favg = fsum / len(flux)
    fuavg = fusum / len(flux_u)
    return favg, fuavg

df = pd.read_csv('psrone_variability.csv')
print(df)
df_154 = df[df["freq (MHz)"] == 154.24]
df_184 = df[df["freq (MHz)"] == 184.96]

fluxes_154 = df_154["flux"]
fluxes_154_u = df_154["flux_error"]
fluxes_184 = df_184["flux"]
fluxes_184_u = df_184["flux_error"]

flux_154, flux_154_u = average_flux(fluxes_154, fluxes_154_u)
flux_184, flux_184_u = average_flux(fluxes_184, fluxes_184_u)

print(f"154 MHz ({len(fluxes_154)}): {flux_154:.2f} +/- {flux_154_u:.2f}")
print(f"184 MHz ({len(fluxes_184)}): {flux_184:.2f} +/- {flux_184_u:.2f}")

freqs     = [154.24, 184.96]
bands     = [30.72, 30.72]
fluxs     = [flux_154, flux_184]
flux_errs = [flux_154_u, flux_184_u]
refs      = ["MWA", "MWA"]

freqs     += [325.0978, 373.9260, 422.7541, 471.5823]
bands     += [50, 50, 50, 50]
fluxs     += [1578.6, 1194.9, 959.3, 591.2]
flux_errs += [223.2, 138.4, 71.4, 98.0]
refs      += ["GMRT_Band_3", "GMRT_Band_3", "GMRT_Band_3", "GMRT_Band_3"]

freqs     += [578.4183, 627.2465, 676.0746, 724.9028]
bands     += [50, 50, 50, 50]
fluxs     += [570.85, 461.35, 550.55, 395.00]
flux_errs += [53.15, 34.35, 31.95, 45.15]
refs      += ["GMRT_Band_4", "GMRT_Band_4", "GMRT_Band_4", "GMRT_Band_4"]

# 700-900, 900-1200, 1200-1600
freqs     += [788.66, 1047.0, 1256.67, 1280.65, 1379.5]
bands     += [200, 300, 128, 128, 400]
fluxs     += [323.16, 1146.3, 90.96, 60.41, 164.02]
flux_errs += [207.92, 454.72, 12.61, 15.35, 46.0]
refs      += ["Parkes", "Parkes", "Parkes", "Parkes", "Parkes"]

find_best_spectral_fit("J0036-1033", freqs, bands, np.array(fluxs), np.array(flux_errs), refs, plot_best=True, alternate_style=True)
