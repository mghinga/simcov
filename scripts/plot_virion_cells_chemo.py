#!/usr/bin/env python3

import time
import os
import sys
import argparse
import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import numpy as np
from sklearn.metrics import r2_score


style.use('fivethirtyeight')
#style.use('qpaper')

plt.rc('axes', titlesize=12)
plt.rc('axes', labelsize=10)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)


argparser = argparse.ArgumentParser()
argparser.add_argument("-f", "--stats-file", required=True, help="Stats file containing simcov results")
#argparser.add_argument("-i", "--indexes", required=True, help="Comma-separated list of column indexes for plotting")
argparser.add_argument("-o", "--output", required=True, help="Output file for images")
argparser.add_argument("-c", "--compare-file", default='', help="File for comparisons")
argparser.add_argument("-p", "--patient-data", default='', help="File with patient data")
argparser.add_argument("-r", "--resolution", type=int, dest='resolution', default=1440, help='Resolution: number of time steps per day') 
#argparser.add_argument("--virus-scale", type=float, dest='virus_scale', default=4e18, help='Factor to scale comparison virus levels')
#argparser.add_argument("--chemo-scale", type=float, dest='chemo_scale', default=5e14, help='Factor to scale comparison chemokine levels')
argparser.add_argument("--cell-scale", type=float, dest='cell_scale', default=1.0, help='Factor to scale comparison cell numbers')
argparser.add_argument("--virus-scale", type=float, dest='virus_scale', default=1.0, help='Factor to scale comparison virus levels')
argparser.add_argument("--chemo-scale", type=float, dest='chemo_scale', default=1.0, help='Factor to scale comparison chemokine levels')
argparser.add_argument("--log", dest='log_scale', action="store_true", help='Use log scale for epicells and tcells')
argparser.add_argument("-d", "--dimensions", required=True, help="Dimensions along one side")
argparser.add_argument("-a", "--animate", action="store_true", help="Animate to see values changing as simcov.stats is update")
options = argparser.parse_args()

#columns = [int(i) for i in options.indexes.split(',')]

fig = plt.figure(figsize=(20, 5))
ax_epicells = fig.add_subplot(1, 3, 2)
ax_tcells = fig.add_subplot(1, 3, 3)
#ax_virus = fig.add_subplot(2, 2, 3)
ax_total_virions = fig.add_subplot(1, 3, 1)

moddate = os.stat(options.stats_file)[8]
unchanged = 0
first = True

def plot_subplot(fname, ax, columns, title, colors, y_label, fig_title, legend_labels, lw=2, alpha=1.0, clear=True, log_scale=False, scale=1.0, xs=None):
    graph_data = open(fname,'r').read()
    lines = graph_data.split('\n')
    get_xs = False
    if xs == None:
        get_xs = True
        xs = []
    ys = []
    labels = []
    fields = lines[0][2:].split('\t')
    for j in range(len(columns)):
        ys.append([])
        labels.append(fields[columns[j]])
    for line in lines[1:]:
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        if get_xs:
            xs.append(float(fields[0]) / options.resolution)
        for j in range(len(columns)):
            if title == 'virions':
                ys[j].append(scale * float(fields[columns[j]]) * options.virus_scale)
            else:
                 ys[j].append(scale * float(fields[columns[j]]))
    if clear:
        ax.clear()
    if log_scale:
        ax.set_yscale('log')
        
    if title == 'patient_data':
        for j in range(len(columns)):
            ax.scatter(xs, ys[j], label='patient_data', lw=0.5, alpha=alpha, color=colors[j], edgecolors='black', zorder=10)
    elif title == 'virions_compare':
        for j in range(len(columns)):
            ax.plot(xs, ys[j], label=labels[j][0:7], lw=lw, alpha=alpha, color=colors[j], zorder=0)
    else:
        for j in range(len(columns)):
            ax.plot(xs, ys[j], label=labels[j][0:7], lw=lw, alpha=alpha, color=colors[j], zorder=5)
    if title == 'epicells':        
        ax.legend(labels=legend_labels, loc='upper left')
    ax.set_xlabel('Time (days)')
    ax.set_title(fig_title)
    if log_scale:
        ax.set_ylim(bottom=1, top=np.power(10, np.ceil(np.log10(np.max(ys))) + 1))
        ax.set_ylabel((y_label + ' in Log10'))
    else:
        ax.set_ylim(bottom=0, top=np.max(ys))
        ax.set_ylabel(y_label)
    xticks = ax.get_xticks()
    if xticks[1] - xticks[0] > 1 and len(xs) > 0:
        ax.set_xticks(range(0, int(max(xs)) + 1, 2))
    plt.tight_layout()
    
    return xs, ys

def animate(i):
    global moddate
    global unchanged
    global first
    new_moddate = os.stat(options.stats_file)[8]
    if new_moddate != moddate or first:
        moddate = new_moddate
        first = False
        
        simcov_xs, simcov_ys = plot_subplot(options.stats_file, ax_epicells, [1, 2, 3, 4], 'epicells',  ['#ffca4f', '#ff3b05', '#197dff', '#000000'], 'Number of Cells in Each State', 'Cells in Each State Over Time', ['Incubating', 'Expressing', 'Apoptotic', 'Dead'], log_scale=options.log_scale)
        #if options.compare_file != '':
            #plot_subplot(options.compare_file, ax_epicells, [2, 3, 5, 4], 'epicells', lw=4, alpha=0.3, clear=False, log_scale=options.log_scale)
            #compare_xs, compare_ys = plot_subplot(options.compare_file, ax_epicells, [1, 2, 3, 4], 'epicells', lw=4, alpha=0.3, clear=False, log_scale=options.log_scale, scale=options.cell_scale)
            
        simcov_xs, simcov_ys = plot_subplot(options.stats_file, ax_tcells, [6], 'tcells', ['#009900'], 'Number of T-Cells', 'T-Cells Over Time', ['T-Cells'], log_scale=options.log_scale)
        #if options.compare_file != '':
            #plot_subplot(options.compare_file, ax_tcells, [6, 7], 'tcells', lw=4, alpha=0.3, clear=False, log_scale=options.log_scale)
            #compare_xs, compare_ys = plot_subplot(options.compare_file, ax_tcells, [6, 5], 'tcells', lw=4, alpha=0.3, clear=False, log_scale=options.log_scale)
            
        simcov_xs, simcov_ys = plot_subplot(options.stats_file, ax_total_virions, [8], 'virions', ['#ff8c00'], 'Number of Virions', 'Virions Over Time', ['Virions'], log_scale=options.log_scale, scale=int(options.dimensions) ** 2)
        #if options.compare_file != '':
        #    compare_xs, compare_ys = plot_subplot(options.compare_file, ax_total_virions, [5], 'virions_compare', lw=4, alpha=0.4, clear=False, log_scale=options.log_scale)
        #if options.patient_data != '':
        #    patient_xs, patient_ys = plot_subplot(options.patient_data, ax_total_virions, [1], 'patient_data', lw=4, alpha=1.0, clear=False, log_scale=options.log_scale)
                
        if options.animate:
            return

    unchanged += 1
    if not options.animate or unchanged > 4:
        moddate = new_moddate
        unchanged = 0
        plt.savefig(options.output + '.pdf', dpi=500, bbox_inches = "tight")
        plt.savefig(options.output + '.png', dpi=500, bbox_inches = "tight")


if options.animate:    
    ani = animation.FuncAnimation(fig, animate, interval=1000)
else:
    animate(0)
#plt.show()

