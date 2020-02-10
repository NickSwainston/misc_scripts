#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.collections as collections

fig, ax = plt.subplots(2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.1},
                       figsize=(20,10))

#fig = plt.figure( figsize=(15,5))
#ax = fig.add_subplot(111)

npointing = 1
secs_plotting = 5

#old timing
#read_time = 0.2 #s
#calc_time = 0.17
#write_time = 0.035

#Galaxy
#retimed with better timing averaging
read_time = 0.516 #s
calc_time = 0.581
write_time = 0.026

#ozstar
read_time = 0.280 #s
cal_time = 0.058
calc_time = 0.095
#calc_time = 0.153 #total calc
write_time = 0.032



#original
read_starts = [0.]
cal_starts = [read_time]
calc_starts = [read_time + cal_time]
write_starts = [read_time + cal_time + calc_time]
for i in range(1,5):
    read_starts.append(write_starts[i-1] + write_time)
    cal_starts.append(read_starts[i] + read_time)
    calc_starts.append(cal_starts[i] + cal_time) 
    write_starts.append(calc_starts[i] + calc_time)

#print(read_starts, calc_starts, write_starts)

linewidth = 1.

for i in range(secs_plotting):
    read  = Rectangle([read_starts[i],  2.],  read_time,  1., facecolor='r',
                      linewidth=linewidth, edgecolor='black')
    ax[0].add_artist(read)
    cal  = Rectangle([cal_starts[i],  1.],  cal_time,  1., facecolor='purple',
                      linewidth=linewidth, edgecolor='black')
    ax[0].add_artist(cal)
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

"""
#multi-pixel
read_starts = [0.]
temp_cal = [read_time]
temp_calc = [read_time + cal_time]
for p in range(1, npointing):
    temp_cal.append(temp_calc[p-1] + calc_time)
    temp_calc.append(temp_cal[p] + cal_time)
cal_starts = [temp_cal]
calc_starts = [temp_calc]
temp_write = [calc_starts[0][-1] + calc_time]
for p in range(1, npointing):
    temp_write.append(temp_write[p-1] + write_time)
write_starts = [temp_write]

for i in range(1, secs_plotting):
    read_starts.append(write_starts[-1][-1] + write_time)
    temp_cal = [read_starts[i] + read_time]
    temp_calc = [temp_cal[0] + cal_time]
    for p in range(1, npointing):
        temp_cal.append(temp_calc[p-1] + calc_time)
        temp_calc.append(temp_cal[p] + cal_time)
    cal_starts.append(temp_cal)
    calc_starts.append(temp_calc) 
    temp_write = [calc_starts[i][-1] + calc_time]
    for p in range(1, npointing):
        temp_write.append(temp_write[p-1] + write_time)
    write_starts.append(temp_write)
    
patches = []
for i in range(secs_plotting):
    read  = Rectangle([read_starts[i],  2.],  read_time, 1., facecolor='r',
                      linewidth=linewidth, edgecolor='black')
    ax[1].add_artist(read)
    for p in range(npointing):
        cal  = Rectangle([cal_starts[i][p],  1.],  cal_time, 1., facecolor='purple',
                          linewidth=linewidth, edgecolor='black')
        ax[1].add_artist(cal)
        calc  = Rectangle([calc_starts[i][p],  1.],  calc_time, 1., facecolor='g',
                          linewidth=linewidth, edgecolor='black')
        ax[1].add_artist(calc)
        write = Rectangle([write_starts[i][p], 0.], write_time, 1., facecolor='b',
                          linewidth=linewidth, edgecolor='black')
        ax[1].add_artist(write)


#temp change of times for async
read_time  *= 5.
cal_time   *= 2.
calc_time  *= 2.
write_time *= 5.

#openMP version:
#read_time = 1.36 #s
#calc_time = 0.22
#write_time = 0.084


#multi pixel async
read_starts = [0.]
temp_cal = [read_time]
temp_calc = [read_time + cal_time]
for p in range(1, npointing):
    temp_cal.append(temp_calc[p-1] + calc_time)
    temp_calc.append(temp_cal[p] + cal_time)
cal_starts = [temp_cal]
calc_starts = [temp_calc]
temp_write = [calc_starts[0][-1] + calc_time]
for p in range(1, npointing):
    temp_write.append(temp_write[p-1] + write_time)
write_starts = [temp_write]

for i in range(1, secs_plotting):
    #async mode
    if i == 1:
        read_starts.append(read_starts[0] + read_time)
    #elif i == 2:
    #    read_starts.append(read_starts[1] + calc_time)
    else:
        read_starts.append(calc_starts[i-2][-1] + calc_time)
    temp_cal = [calc_starts[i-1][-1] + calc_time]
    temp_calc = [temp_cal[0] + cal_time]
    for p in range(1, npointing):
        temp_cal.append(temp_calc[p-1] + calc_time)
        temp_calc.append(temp_cal[p] + cal_time)
    cal_starts.append(temp_cal)
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
    ax[2].add_artist(read)
    for p in range(npointing):
        cal  = Rectangle([cal_starts[i][p],  1.],  cal_time, 1., facecolor='purple',
                          linewidth=linewidth, edgecolor='black')
        ax[2].add_artist(cal)
        calc  = Rectangle([calc_starts[i][p],  1.],  calc_time, 1., facecolor='g',
                          linewidth=linewidth, edgecolor='black')
        ax[2].add_artist(calc)
        write = Rectangle([write_starts[i][p], 0.], write_time, 1., facecolor='b',
                          linewidth=linewidth, edgecolor='black')
        ax[2].add_artist(write)
#p = collections.PatchCollection(patches)
#ax.add_collection(p)

max_xaxis = write_starts[-1][-1] + write_time
"""


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

#Ozstar
#serial and cal once 15p:
read_time  = 0.266 #s
calc_time  = 0.044
write_time = 0.013

#Ozstar
#serial and cal once 1p:
read_time  = 0.266 #s
cal_time   = 0.070
#calc_time  = 0.085 Calculated using calc_time_total * .55
calc_time  = 0.039 # calculated from the (15p_calc_time - 0.07)/15

#calc_time = 0.155 #total
write_time = 0.013




"""
#remove unnessary calcs
read_starts = [0.]
temp_cal = [read_time]
temp_calc = [read_time + cal_time]
for p in range(1, npointing):
    temp_cal.append(temp_calc[p-1] + calc_time)
    temp_calc.append(temp_cal[p] + cal_time)
cal_starts = [temp_cal]
calc_starts = [temp_calc]
temp_write = [calc_starts[0][-1] + calc_time]
for p in range(1, npointing):
    temp_write.append(temp_write[p-1] + write_time)
write_starts = [temp_write]

for i in range(1, secs_plotting):
    read_starts.append(write_starts[-1][-1] + write_time)
    temp_cal = [read_starts[i] + read_time]
    temp_calc = [temp_cal[0] + cal_time]
    for p in range(1, npointing):
        temp_cal.append(temp_calc[p-1] + calc_time)
        temp_calc.append(temp_cal[p] + cal_time)
    cal_starts.append(temp_cal)
    calc_starts.append(temp_calc) 
    temp_write = [calc_starts[i][-1] + calc_time]
    for p in range(1, npointing):
        temp_write.append(temp_write[p-1] + write_time)
    write_starts.append(temp_write)
    
patches = []
for i in range(secs_plotting):
    read  = Rectangle([read_starts[i],  2.],  read_time, 1., facecolor='r',
                      linewidth=linewidth, edgecolor='black')
    ax[3].add_artist(read)
    for p in range(npointing):
        cal  = Rectangle([cal_starts[i][p],  1.],  cal_time, 1., facecolor='purple',
                          linewidth=linewidth, edgecolor='black')
        ax[3].add_artist(cal)
        calc  = Rectangle([calc_starts[i][p],  1.],  calc_time, 1., facecolor='g',
                          linewidth=linewidth, edgecolor='black')
        ax[3].add_artist(calc)
        write = Rectangle([write_starts[i][p], 0.], write_time, 1., facecolor='b',
                          linewidth=linewidth, edgecolor='black')
        ax[3].add_artist(write)

"""
#multi-pixel
read_starts = [0.]
cal_starts = [read_time]
temp_calc = [read_time + cal_time]
for p in range(1, npointing):
    temp_calc.append(temp_calc[p-1] + calc_time)
calc_starts = [temp_calc]
temp_write = [read_time + cal_time + calc_time*npointing]
for p in range(1, npointing):
    temp_write.append(temp_write[p-1] + write_time)
write_starts = [temp_write]
for i in range(1, secs_plotting):
    read_starts.append(write_starts[-1][-1] + write_time)
    cal_starts.append(read_starts[i] + read_time)
    temp_calc = [cal_starts[i] + cal_time]
    for p in range(1, npointing):
        temp_calc.append(temp_calc[p-1] + calc_time)
    calc_starts.append(temp_calc) 
    temp_write = [calc_starts[i][-1] + calc_time]
    for p in range(1, npointing):
        temp_write.append(temp_write[p-1] + write_time)
    write_starts.append(temp_write)

max_xaxis = write_starts[-1][-1] + write_time

#print(read_starts, calc_starts, write_starts)

patches = []
for i in range(secs_plotting):
    read  = Rectangle([read_starts[i],  2.],  read_time, 1., facecolor='r',
                      linewidth=linewidth, edgecolor='black')
    ax[1].add_artist(read)
    cal  = Rectangle([cal_starts[i],  1.],  cal_time, 1., facecolor='purple',
                      linewidth=linewidth, edgecolor='black')
    ax[1].add_artist(cal)
    for p in range(npointing):
        calc  = Rectangle([calc_starts[i][p],  1.],  calc_time, 1., facecolor='g',
                          linewidth=linewidth, edgecolor='black')
        ax[1].add_artist(calc)
        write = Rectangle([write_starts[i][p], 0.], write_time, 1., facecolor='b',
                          linewidth=linewidth, edgecolor='black')
        ax[1].add_artist(write)
#p = collections.PatchCollection(patches)
#ax.add_collection(p)


#for i, label in enumerate(('a)', 'b)', 'c)', 'd)', 'e)')):
#    ax[i].text(-0.1, 0.5, label, transform=ax[i].transAxes,
#      fontsize=25, va='center', ha='right')


plt.axis([0, max_xaxis, 0, 3])
plt.yticks([0.5, 1.5, 2.5], ['Write', 'Calc', 'Read'])
for axi in range(2):
    for tick in ax[axi].yaxis.get_major_ticks():
        tick.label.set_fontsize(30)
plt.xticks(fontsize=30)
plt.xlabel('Processing time (s)', fontsize=30)
plt.savefig('multi-pixel_workflow.eps', bbox_inches='tight', dpi=1000)


