#script for DM plot

import numpy as np
from numpy.linalg import lstsq
from numpy import genfromtxt 
from numpy import math
from scipy import optimize

# MatPlotlib
import matplotlib.pyplot as plt
from matplotlib import pylab

my_data = genfromtxt('DM_final_plot.csv', delimiter=',')
x = my_data.T[0]
y = my_data.T[2]
z = my_data.T[4]
t = my_data.T[5]


x1 = my_data.T[0][0:1]
y1 = my_data.T[2][0:1]
z1 = my_data.T[4][0:1]
t1 = my_data.T[5][0:1]

x2 = my_data.T[0][1:2]
y2 = my_data.T[2][1:2]
z2 = my_data.T[4][1:2]
t2 = my_data.T[5][1:2]

x3 = my_data.T[0][2:3]
y3 = my_data.T[2][2:3]
z3 = my_data.T[4][2:3]
t3 = my_data.T[5][2:3]

x4 = my_data.T[0][3:4]
y4 = my_data.T[2][3:4]
z4 = my_data.T[4][3:4]
t4 = my_data.T[5][3:4]

x5 = my_data.T[0][4:5]
y5 = my_data.T[2][4:5]
z5 = my_data.T[4][4:5]
t5 = my_data.T[5][4:5]

x6 = my_data.T[0][5:6]
y6 = my_data.T[2][5:6]
z6 = my_data.T[4][5:6]
t6 = my_data.T[5][5:6]

x7 = my_data.T[0][6:7]
y7 = my_data.T[2][6:7]
z7 = my_data.T[4][6:7]
t7 = my_data.T[5][6:7]

x8 = my_data.T[0][7:8]
y8 = my_data.T[2][7:8]
z8 = my_data.T[4][7:8]
t8 = my_data.T[5][7:8]

logx = np.log10(x)
logy = np.log10(y)
logz = np.log10(z)


def power_law(x, a, b, c):
    return a*np.power(x, b) + c
    
# Fit the power-law to data
params, cov = optimize.curve_fit(power_law, x, y, p0=np.asarray([-0.00001, 3.68, 11.41]), sigma=z, maxfev=5000, bounds=(-np.inf, np.inf))
print(params)

#plt.errorbar(x, abs(y), yerr=z, xerr=t, fmt='.', ecolor='white', color='white', capsize=0.1)
plt.errorbar(x1, abs(y1), yerr=z1, xerr=t1, fmt='.', markersize=5, color='red', ecolor='black', capsize=4, label='MWA')
plt.errorbar(x2, abs(y2), yerr=z2, xerr=t2, fmt='.', markersize=5, color='blue', ecolor='black', capsize=4, label='GMRT')
plt.errorbar(x3, abs(y3), yerr=z3, xerr=t3, fmt='.', markersize=5, color='yellow', ecolor='black', capsize=4, label='Parkes B1')
plt.errorbar(x4, abs(y4), yerr=z4, xerr=t4, fmt='.', markersize=5, color='magenta', ecolor='black', capsize=4, label='Parkes B2')
plt.errorbar(x5, abs(y5), yerr=z5, xerr=t5, fmt='.', markersize=5, color='cyan', ecolor='black', capsize=4, label='Parkes B3')
plt.errorbar(x6, abs(y6), yerr=z6, xerr=t6, fmt='.', markersize=5, color='yellow', ecolor='black', capsize=4)
plt.errorbar(x7, abs(y7), yerr=z7, xerr=t7, fmt='.', markersize=5, color='magenta', ecolor='black', capsize=4)
plt.errorbar(x8, abs(y8), yerr=z8, xerr=t8, fmt='.', markersize=5, color='cyan', ecolor='black', capsize=4)

#plt.ylim(-0.005, 0.02)
plt.yscale('log')
plt.xscale('log')
x = np.arange(min(x), max(x), 1)
plt.plot(x, power_law(x, params[0], params[1], params[2]), 'k', linestyle='--', label='label name')
#plt.plot(x, power_law(x, -300, 10, 11.41), 'k', linestyle='--', label='label name')

plt.title('07-09 Nov 2019')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Measured DM (pc cm-3)")
plt.tight_layout(pad=0.5)    
   

plt.savefig('try.png')
#plt.show()
