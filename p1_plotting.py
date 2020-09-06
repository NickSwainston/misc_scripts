import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x, a, b):
    return a * x + b

# The data: gps_time = [P] and absorbance, period
#gps_time = np.array([1137236608, 1150234552, 1164110416, 1182630616, 1194350120, 1225462936, 1222697776])
#period   = np.array([0.9000079,  0.9000096,  0.9000090,  0.9000382,  0.8999792,  0.900019,   0.90000818])
#period_u = np.array([0.0000029,  0.0000012,  0.0000014,  0.0000039,  0.0000016,  0.0000045,  0.00000401])

#archive obs
gps_time = np.array([1222697776, 1137236608, 1194350120, 1225462936, 1182630616, 1220886016, 1133775752, 1220886016, 1255444104])
period   = np.array([900.0046, 900.011, 900.012, 900.01, 900.006, 900.009, 900.011, 900.004, 900.007]) / 1000
period_u = np.array([0.004, 0.005, 0.006, 0.006, 0.008, 0.008, 0.013, 0.008, 0.002]) / 1000
freq = 1./ period
freq_u = 1./period_u

#DDT obs
gps_time = np.array([1274143152, 1274573984, 1275085816, 1275092416, 1275094456, 1275172216, 1275177136, 1275178816, 1275258616, 1275431416, 1275863416, 1275866536, 1276725752, 1278106408])
period   = np.array([900.0010936742, 900.0207476841, 900.0101253760, 900.0094419584, 900.0062279132, 900.0042987105, 899.9992088839, 900.0176455286, 900.0196082668, 900.0058411118, 900.0045834992, 900.0062400582, 900.0111874337, 900.0125086498])/1000
period_u = np.array([0.0079617690, 0.0166588190, 0.0082400282, 0.0082397255, 0.0082400535, 0.0082399215, 0.0080716686, 0.0082402626, 0.0082398148, 0.0053934849, 0.0082399267, 0.0080719841, 0.0053445405, 0.0053933117]) /1000

popt, pcov = curve_fit(line, gps_time, period, sigma=period_u)
fopt, fcov = curve_fit(line, gps_time, 1/period, sigma=period_u/period**2)
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

