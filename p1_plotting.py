import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x, a, b):
    return a * x + b

# The data: gps_time = [P] and absorbance, period
gps_time = np.array([1137236608, 1150234552, 1164110416, 1182630616, 1194350120, 1225462936, 1222697776])
period   = np.array([0.9000079,  0.9000096,  0.9000090,  0.9000382,  0.8999792,  0.900019,   0.90000818])
period_u = np.array([0.0000029,  0.0000012,  0.0000014,  0.0000039,  0.0000016,  0.0000045,  0.00000401])
freq = 1./ period
freq_u = 1./period_u

popt, pcov = curve_fit(line, gps_time, period, sigma=period_u)
print("P1 =", popt[0], "+/-", pcov[0,0]**0.5)
fopt, fcov = curve_fit(line, gps_time, freq, sigma=freq_u)
print("F1 =", fopt[0], "+/-", fcov[0,0]**0.5)

fig, ax = plt.subplots()
plt.title("MWA pulsar candidate DM 23.039 P 900 ms\nP1 = {0} +/- {1} s/s".format(popt[0], pcov[0,0]**0.5))
plt.errorbar(gps_time, period, yerr=period_u, fmt=".")
plt.plot(gps_time, line(gps_time, popt[0], popt[1]), 'r-')
ax.ticklabel_format(useOffset=False)
plt.xlabel('GPS time (s)')
plt.ylabel('Period (s)')
#plt.show()
plt.savefig('p1_fit.png', bbox_inches='tight')

plt.clr()

plt.title("MWA pulsar candidate DM 23.039 P 900 ms\nF1 = {0} +/- {1} s/s".format(fopt[0], fov[0,0]**0.5))
plt.errorbar(gps_time, freq, yerr=freq_u, fmt=".")
plt.plot(gps_time, line(gps_time, fopt[0], fopt[1]), 'r-')
ax.ticklabel_format(useOffset=False)
plt.xlabel('GPS time (s)')
plt.ylabel('Freq (s^-1)')
#plt.show()
plt.savefig('f1_fit.png', bbox_inches='tight')

