#!/usr/bin/env python

import numpy as np
from numpy.linalg import lstsq
from numpy import genfromtxt 
from numpy import math
from scipy import optimize
from scipy.special import huber
from scipy.optimize import least_squares
from argparse import Namespace

# MatPlotlib
import matplotlib.pyplot as plt
from matplotlib import pylab
from matplotlib.ticker import *

from iminuit import Minuit
from iminuit.cost import LeastSquares


my_data = genfromtxt('plot_candidate_spectral_index.csv', delimiter=',')
freq = my_data.T[0] # MHz
flux = my_data.T[1] / 1000 #mJy
flux_err = my_data.T[2] / 1000 #mJy


logx = np.log10(freq)
logy = np.log10(flux)
logz = flux_err/flux#np.log10(flux_err)

def gen_data(t, a, b, c, noise=0, n_outliers=0, random_state=0):

    y = a + b * np.exp(t * c)


    rnd = np.random.RandomState(random_state)

    error = noise * rnd.randn(t.size)

    outliers = rnd.randint(0, t.size, n_outliers)

    error[outliers] *= 10


    return y + error

# our line model, unicode parameter names are supported :) 
def line_func(x, α, β):
    return β + x * α
def fun(x, t, y):
    return x[1] + x[0] * t - y


def power_law_func(x, α, β):
    return β * x ** α

# Huber loss fucntion
def huber_loss(t, k=1.345):
    """
    if abs(t) < k:
        return 0.5 * t**2
    else:
        return k * abs(t) - 0.5 * k**2
    """
    if abs(t) < k:
        return t
    else:
        return 2 * k * abs(t)**0.5 - k**2
# vectroise it
from numpy import vectorize
v_huber_loss = vectorize(huber_loss)

fits = {}
# Perform fit using different methods
# Curve fit method without errors
params, params_covariance = optimize.curve_fit(line_func, logx, logy, p0=[2,2])
perr = np.sqrt(np.diag(params_covariance))
print('Curve fit method without errors')
print('alpha = {:6.3f}+/- {:6.3f}'.format(params[0], perr[0]))
fits['Curve fit method without errors'] = [params, perr]

# Curve fit method with errors
params, params_covariance = optimize.curve_fit(line_func, logx, logy, sigma=logz, p0=[2,2])
perr = np.sqrt(np.diag(params_covariance))
print('Curve fit method with errors')
print('alpha = {:6.3f}+/- {:6.3f}'.format(params[0], perr[0]))
fits['Curve fit method with errors'] = [params, perr]

# Curve fit method with errors in linear space
params, params_covariance = optimize.curve_fit(power_law_func, freq, flux, sigma=flux_err, p0=[2,2])
perr = np.sqrt(np.diag(params_covariance))
print('Curve fit method with errors in linear space')
print('alpha = {:6.3f}+/- {:6.3f}'.format(params[0], perr[0]))
fits['Curve fit method with errors in linear space'] = [params, perr]

# Scipy least squares method
x0 = np.array([1.0, 1.0])
res_lsq = least_squares(fun, x0, args=(logx, logy), verbose=0)
print('Scipy least squares method')
print('alpha = {:6.3f}'.format(res_lsq['x'][0]))
fits['Scipy least squares method'] = res_lsq['x']

# Scipy huber loss method
x0 = np.array([1.0, 1.0])
res_lsq = least_squares(fun, x0, args=(logx, logy), loss="huber", verbose=0)
print('Scipy huber loss method')
print('alpha = {:6.3f}'.format(res_lsq['x'][0]))
fits['Scipy huber loss method'] = res_lsq['x']


# Minuit least squares method
least_squares = LeastSquares(logx, logy, logz, line_func)
line_func.errordef = Minuit.LEAST_SQUARES
m = Minuit(least_squares, α=-1.5, β=4.5)
m.migrad()  # run optimiser
params = list(m.values)
m.hesse()   # run covariance estimator
perr = list(m.errors)
print('Minuit least squares method')
print('alpha = {:6.3f}+/- {:6.3f}'.format(params[0], perr[0]))
fits['Minuit least squares method'] = [params, perr]

# Minuit least squares method in linear space
least_squares = LeastSquares(freq, flux, flux_err,power_law_func)
power_law_func.errordef = Minuit.LEAST_SQUARES
m = Minuit(least_squares, α=-1.5, β=4.5)
m.migrad()  # run optimiser
params = list(m.values)
m.hesse()   # run covariance estimator
perr = list(m.errors)
print('Minuit least squares method in linear space')
print('alpha = {:6.3f}+/- {:6.3f}'.format(params[0], perr[0]))
fits['Minuit least squares method in linear space'] = [params, perr]

# Soft l1 loss method
least_squares = LeastSquares(logx, logy, logz, line_func, loss="soft_l1")
line_func.errordef = Minuit.LEAST_SQUARES
m = Minuit(least_squares, α=-1.5, β=4.5)
m.migrad()  # run optimiser
params = list(m.values)
m.hesse()   # run covariance estimator
perr = list(m.errors)
print('soft_l1 loss method')
print('alpha = {:6.3f}+/- {:6.3f}'.format(params[0], perr[0]))
fits['soft_l1 loss method'] = [params, perr]

# Huber loss method
least_squares = LeastSquares(logx, logy, logz, line_func, loss=v_huber_loss)
line_func.errordef = Minuit.LEAST_SQUARES
m = Minuit(least_squares, α=-1.5, β=4.5)
m.migrad()  # run optimiser
params = list(m.values)
m.hesse()   # run covariance estimator
perr = list(m.errors)
print('Huber loss method')
print('alpha = {:6.3f}+/- {:6.3f}'.format(params[0], perr[0]))
fits['Huber loss method'] = [params, perr]

#print(m.errors)
"""
print(' ')
print('parameters (a & b)')
print(params)
#print(' ')
#print('covariance (a & b)')
#print(params_covariance)
print(' ')
print('print error (a & b)')
print(perr)
print(' ')
"""

# Split into bands
G3_f = []
G3_d = []
G3_e = []
G4_f = []
G4_d = []
G4_e = []
P_f = []
P_d = []
P_e = []
for i in range(len(freq)):
    if freq[i] < 500:
        # GMRT band 3
        G3_f.append(freq[i]/1000)
        G3_d.append(flux[i])
        G3_e.append(flux_err[i])
    elif freq[i] > 750:
        # Parkes
        P_f.append(freq[i]/1000)
        P_d.append(flux[i])
        P_e.append(flux_err[i])
    else:
        # GMRT band 4
        G4_f.append(freq[i]/1000)
        G4_d.append(flux[i])
        G4_e.append(flux_err[i])

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)

markersize=5
# Plot each sub band in a different colour
plt.errorbar(0.150, 12.5, yerr=2, uplims=True, fmt='.', ecolor='r',
             label='MWA 150 MHz', color='r', capsize=3,  markersize=0.01)
plt.errorbar(G3_f, G3_d, yerr=G3_e, fmt='^', ecolor='g',
             label='GMRT Band 3',color='g', capsize=3,  markersize=markersize)
plt.errorbar(G4_f, G4_d, yerr=G4_e, fmt='D', ecolor='b',
             label='GMRT Band 4',color='b', capsize=3,  markersize=markersize)
plt.errorbar(P_f, P_d, yerr=P_e, fmt='o', ecolor='m',
             label='Parkes UWL',color='m', capsize=3,  markersize=markersize)

# Plot model fit
model_key = 'Curve fit method with errors in linear space'
plt.plot(np.array(freq)/1000, power_law_func(freq, fits[model_key][0][0], fits[model_key][0][1]), 'k',
         linestyle='--', label=r'$\alpha$ = '+u"\u2212"+'{:2.1f}'.format(abs(fits[model_key][0][0]))+u"\u00B1"+'{:2.1f}'.format(abs(fits[model_key][1][0])))
"""
for model_key in fits.keys():
    #print(fits[model_key][0])
    plt.plot(np.array(freq)/1000, 10**(line_func(logx, fits[model_key][0], fits[model_key][1])),
             linestyle='--', label=r'{} $\alpha$ = '.format(model_key)+u"\u2212"+'{:3.2f}'.format(abs(float(fits[model_key][0]))))
"""
#print(np.array([150, freq[0]])/1000)
#print(10**(line_func(np.array([150., float(freq[0])]), params[0], params[1])))
#print(10**(line_func(logx, params[0], params[1])))
log_x_mwa = np.log10([150, freq[0]])
plt.plot(np.array([150, freq[0]])/1000, power_law_func([150, freq[0]], fits[model_key][0][0], fits[model_key][0][1]),
         'k', linestyle='dotted')
"""
plt.errorbar(logx, logy, yerr=logz, fmt='o', ecolor='m',
             label='Parkes UWL',color='m', capsize=3,  markersize=3)
plt.plot(logx, line_func(logx, params[0], params[1]), 'k',
         linestyle='--', label=r'$\alpha$ = '+u"\u2212"+'{:3.2f}'.format(abs(params[0])))

"""
print("Outlier t value {}".format((line_func(np.log10(1047.0), params[0], params[1]) - np.log10(1146.3)) / (454.72/1146.3)))
plt.yscale('log')
plt.xscale('log')
#ax.semilogy(range(10))

ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_major_formatter(ScalarFormatter())
x_minor_range = np.array(range(200,1000,200))/1000
#print(x_minor_range)
ax.xaxis.set_minor_locator(FixedLocator(x_minor_range))
ax.xaxis.set_minor_formatter(FixedFormatter(x_minor_range))

ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))


ax.tick_params(top=True, right=True, which='both', labelsize=12)
#ax.tick_params(axis="x", labelrotation=-45, which='both') 
plt.legend()
plt.xlabel("Frequency (GHz)", fontsize=12)
plt.ylabel("Flux density (mJy)", fontsize=12)
plt.tight_layout()
plt.savefig('plot_candidate_spectral_index.png', dpi=200)
plt.savefig('plot_candidate_spectral_index.eps')
