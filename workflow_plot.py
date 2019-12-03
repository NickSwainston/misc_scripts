#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.collections as collections


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


