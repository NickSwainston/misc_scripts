#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.collections as collections
"""
n_pointing = np.arange(1,31)

#Galaxy data---------------------------------------
#original
G_orig    = np.full(len(n_pointing), 1.636) * n_pointing * 24.
G_orig_er = np.full(len(n_pointing), 0.048) * 24.

#Multi-pixel beamformer
#TODO updata with pointing range
G_mpb    = np.array([2.831392019, 1.600111409, 1.256015202, 1.143495984, 1.132492969,
                     1.137000819, 1.146377474, 1.1315021,   1.13003097,  1.124950794,
                     1.12363004,  1.122008514, 1.117132084, 1.120854555, 1.115692496,
                     1.137000819, 1.146377474, 1.1315021,   1.13003097,  1.124950794,
                     1.12363004,  1.122008514, 1.117132084, 1.120854555, 1.115692496,
                     1.117132084, 1.120854555, 1.115692496, 1.137000819, 1.146377474])\
                     * n_pointing * 24.
G_mpb_er = np.array([0.08228635476, 0.02443385882, 0.01424811249, 0.004902127501, 0.002967214197,
                     0.01245005442, 0.02087704667, 0.00381140615, 0.006075820599, 0.004261909977,
                     0.00337042659, 0.00186423623, 0.00224264162, 0.003157883716, 0.003126045247,
                     0.01245005442, 0.02087704667, 0.00381140615, 0.006075820599, 0.004261909977,
                     0.00337042659, 0.00186423623, 0.00224264162, 0.003157883716, 0.003126045247,
                     0.00381140615, 0.006075820599, 0.0042619099, 0.00337042659, 0.00186423623])\
                     * 24.

#Ozstar data---------------------------------------
#original
O_orig    = np.full(len(n_pointing), 0.736) * n_pointing * 24.
O_orig_er = np.full(len(n_pointing), 0.084) * 24.

#Multi-pixel beamformer
#TODO updata with pointing range
O_mpb    = np.array([3.139471429,  1.516942857, 0.9796420635, 0.6210142857, 0.6109428571,
                     0.4121063492, 0.398287415, 0.3594440476, 0.2960134921, 0.2475538095,
                     0.2379311688, 0.207063293, 0.203425641,  0.2021408163, 0.2136314286,
                     0.2115669643, 0.207242577, 0.2067165344, 0.2238288221, 0.2023686905,
                     0.2168582766, 0.205066774, 0.2237901656, 0.2181756944, 0.2171654286,
                     0.2170190476, 0.215439418, 0.2140656463, 0.2138665025, 0.2168870635])\
                     * n_pointing * 24.
O_mpb_er = np.array([0.3031622606, 0.1867256804, 0.1189938961, 0.06173374331, 0.07661494931,
                     0.0551306001, 0.0589499940, 0.0490239472, 0.03659994481, 0.01882398827,
                     0.0188806783, 0.0073943070, 0.0113471232, 0.00859965614, 0.00997227971,
                     0.0056918661, 0.0070319188, 0.0114938482, 0.01259777255, 0.00665450730,
                     0.0124219095, 0.0080223915, 0.0113834284, 0.00772366071, 0.00841511883,
                     0.0071728955, 0.0098200043, 0.0068575750, 0.00829136570, 0.00819338188])\
                     * n_pointing * 24.


plt.errorbar(n_pointing, G_orig, yerr=G_orig_er, label='Galaxy original beamformer')
plt.errorbar(n_pointing, G_mpb,  yerr=G_mpb_er,  label='Galaxy multi-pixel beamformer')
plt.errorbar(n_pointing, O_orig, yerr=O_orig_er, label='Ozstar original beamformer')
plt.errorbar(n_pointing, O_mpb,  yerr=O_mpb_er,  label='Ozstar multi-pixel beamformer')

plt.xlabel('Number of tied-array beams')
plt.ylabel('Processing time (s) per second of data')
plt.legend(loc='upper left')
#plt.show()
plt.savefig("benchmarking_total.eps")

plt.clf()

plt.errorbar(n_pointing, G_orig / n_pointing, yerr=G_orig_er, label='Galaxy original beamformer')
plt.errorbar(n_pointing, G_mpb / n_pointing,  yerr=G_mpb_er,  label='Galaxy multi-pixel beamformer')
plt.errorbar(n_pointing, O_orig / n_pointing, yerr=O_orig_er, label='Ozstar original beamformer')
plt.errorbar(n_pointing, O_mpb / n_pointing,  yerr=O_mpb_er / n_pointing,  label='Ozstar multi-pixel beamformer')

plt.xlabel('Number of tied-array beams')
plt.ylabel('Processing time (s) per second of data')
plt.legend(loc='upper left')
#plt.show()
plt.savefig("benchmarking_per_pointing.eps")
"""
# Make workflow plots---------------------------------------------------------------------

plt.clf()
fig, ax = plt.subplots(2, sharex=True, sharey=True, gridspec_kw={'hspace': 0}, figsize=(15,10))

#fig = plt.figure( figsize=(15,5))
#ax = fig.add_subplot(111)

npointing = 1
secs_plotting = 5

#old timing
#read_time = 0.2 #s
#calc_time = 0.17
#write_time = 0.035

#retimed with better timing averaging
read_time = 0.516 #s
calc_time = 0.581
write_time = 0.026



#original
read_starts = [0.]
calc_starts = [read_time]
write_starts = [npointing*calc_time + read_time]
for i in range(1,5):
    read_starts.append(write_starts[i-1] + write_time)
    calc_starts.append(read_starts[i] + read_time) 
    write_starts.append(calc_starts[i] + calc_time)

#print(read_starts, calc_starts, write_starts)

linewidth = 2.

for i in range(secs_plotting):
    read  = Rectangle([read_starts[i],  2.],  read_time,  1., facecolor='r',
                      linewidth=linewidth, edgecolor='black')
    ax[0].add_artist(read)
    calc  = Rectangle([calc_starts[i],  1.],  calc_time,  1., facecolor='g',
                      linewidth=linewidth, edgecolor='black')
    ax[0].add_artist(calc)
    write = Rectangle([write_starts[i], 0.], write_time, 1., facecolor='b',
                      linewidth=linewidth, edgecolor='black')
    ax[0].add_artist(write)
#plt.axis([0, write_starts[-1]+write_time, 0, 3])
#plt.yticks([0.5, 1.5, 2.5], ['Write', 'Calc', 'Read'], fontsize=20)
#plt.xlabel('Processing time (s)')
#plt.savefig('original_workflow.eps')
#plt.yticks([0.5, 1.5, 2.5], ['Write', 'Calc', 'Read'], fontsize=30)

#plt.clf()
#plt.close()
#fig = plt.figure(figsize=(15,5))
#ax = fig.add_subplot(111)

npointing = 5

#openMP version:
#read_time = 1.36 #s
#calc_time = 0.22
#write_time = 0.084

#serial upgrade 1p:
read_time  = 0.542 #s
calc_time  = 0.510
write_time = 0.009

#serial upgrade 15p:
read_time  = 0.471 #s
calc_time  = 0.368
write_time = 0.012

#serial and cal once 1p:
read_time  = 0.499 #s
calc_time  = 0.363
write_time = 0.008

#serial and cal once 15p:
read_time  = 0.516 #s
calc_time  = 0.275
write_time = 0.016


#multi-pixel
read_starts = [0.]
temp_calc = [read_time]
for p in range(1, npointing):
    temp_calc.append(temp_calc[p-1] + calc_time)
calc_starts = [temp_calc]
temp_write = [read_time + calc_time*npointing]
for p in range(1, npointing):
    temp_write.append(temp_write[p-1] + write_time)
write_starts = [temp_write]
for i in range(1, secs_plotting):
    """
    #async mode
    if i == 1:
        read_starts.append(read_starts[0] + read_time)
    #elif i == 2:
    #    read_starts.append(read_starts[1] + calc_time)
    else:
        read_starts.append(calc_starts[i-2][0] + npointing*calc_time)
    temp_calc = [calc_starts[i-1][-1] + calc_time]
    for p in range(1, npointing):
        temp_calc.append(temp_calc[p-1] + calc_time)
    calc_starts.append(temp_calc) 
    temp_write = [calc_starts[i][-1] + calc_time]
    for p in range(1, npointing):
        temp_write.append(temp_write[p-1] + write_time)
    write_starts.append(temp_write)
    """
    read_starts.append(calc_starts[-1][-1] + calc_time)
    temp_calc = [read_starts[i] + read_time]
    for p in range(1, npointing):
        temp_calc.append(temp_calc[p-1] + calc_time)
    calc_starts.append(temp_calc) 
    temp_write = [calc_starts[i][-1] + calc_time]
    for p in range(1, npointing):
        temp_write.append(temp_write[p-1] + write_time)
    write_starts.append(temp_write)


#print(read_starts, calc_starts, write_starts)

patches = []
for i in range(secs_plotting):
    read  = Rectangle([read_starts[i],  2.],  read_time, 1., facecolor='r',
                      linewidth=linewidth, edgecolor='black')
    ax[1].add_artist(read)
    for p in range(npointing):
        calc  = Rectangle([calc_starts[i][p],  1.],  calc_time, 1., facecolor='g',
                          linewidth=linewidth, edgecolor='black')
        ax[1].add_artist(calc)
    for p in range(npointing):
        write = Rectangle([write_starts[i][p], 0.], write_time, 1., facecolor='b',
                          linewidth=linewidth, edgecolor='black')
        ax[1].add_artist(write)
#p = collections.PatchCollection(patches)
#ax.add_collection(p)
print(write_starts[0][-1], write_starts[1][-1], write_starts[2][-1], write_starts[3][-1], write_starts[4][-1])
plt.axis([0, write_starts[-1][-1] + write_time, 0, 3])
plt.yticks([0.5, 1.5, 2.5], ['Write', 'Calc', 'Read'])
for tick in ax[0].yaxis.get_major_ticks():
    tick.label.set_fontsize(30)
for tick in ax[1].yaxis.get_major_ticks():
    tick.label.set_fontsize(30)
plt.xticks(fontsize=30)
plt.xlabel('Processing time (s)', fontsize=30)
plt.savefig('multi-pixel_workflow.eps')


