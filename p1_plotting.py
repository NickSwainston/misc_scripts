import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x, a, b):
    return a * x + b

# The data: gps_time = [P] and absorbance, period
#gps_time = np.array([1137236608, 1150234552, 1164110416, 1182630616, 1194350120, 1225462936, 1222697776])
#period   = np.array([0.9000079,  0.9000096,  0.9000090,  0.9000382,  0.8999792,  0.900019,   0.90000818])
#period_u = np.array([0.0000029,  0.0000012,  0.0000014,  0.0000039,  0.0000016,  0.0000045,  0.00000401])

gps_time = np.array([1222697776, 1137236608, 1194350120, 1225462936, 1182630616, 1220886016, 1133775752, 1220886016, 1255444104])
period   = np.array([900.0046, 900.011, 900.012, 900.01, 900.006, 900.009, 900.011, 900.004, 900.007]) / 1000
period_u = np.array([0.004, 0.005, 0.006, 0.006, 0.008, 0.008, 0.013, 0.008, 0.002]) / 1000

popt, pcov = curve_fit(line, gps_time, period, sigma=period_u)
print("P1 =", popt[0], "+/-", pcov[0,0]**0.5)

fig, ax = plt.subplots()
plt.title("MWA pulsar candidate DM 23.039 P 900 ms\nP1 = {0} +/- {1} s/s".format(popt[0], pcov[0,0]**0.5))
plt.errorbar(gps_time, period, yerr=period_u, fmt=".")
plt.plot(gps_time, line(gps_time, popt[0], popt[1]), 'r-')
ax.ticklabel_format(useOffset=False)
plt.xlabel('GPS time (s)')
plt.ylabel('Period (s)')
#plt.show()
plt.savefig('p1_fit.png', bbox_inches='tight')
