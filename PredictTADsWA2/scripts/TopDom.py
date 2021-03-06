#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
import signal_func
import pandas as pd
import argparse


def TopDom(file, window, chrom, res, peak_find_funk=None,
           smooth_func=False, filter_func=False, bin_signal=False):
    matrix = np.genfromtxt(file, delimiter='\t')
    s = matrix.shape
    assert s[0] == s[1], 'Unknown Type of matrix file'
    n_bins = s[0]
    bin_sig = binSignal(matrix, window, n_bins)
    no_gap_region = find_not_gap(matrix, window)
    all_peaks = []
    df = []
    for i in no_gap_region:
        print('process:', i[0], i[1])
        ana_reg = bin_sig[i[0]:i[1]]
        if smooth_func is not False:
            ana_reg = smooth_func(ana_reg)
        peaks = peak_find_funk(ana_reg)
        if filter_func is not False:
            if filter_func != signal_func.statFilter:
                peaks = filter_func(ana_reg, peaks)
        peaks = peaks + i[0]
        all_peaks.extend(peaks)
    if filter_func == signal_func.statFilter:
        all_peaks, p_values = signal_func.statFilter(matrix,
                                                     window, no_gap_region,
                                                     all_peaks)
    all_peaks = list(all_peaks)
    all_peaks.extend([no_gap_region[i][1] for i in range(len(no_gap_region))])
    df = pd.DataFrame(all_peaks)
    df['type'] = 'domain'
    df2 = [no_gap_region[i][0]-1 for i in range(len(no_gap_region))]
    df2 = pd.DataFrame(df2)
    df2['type'] = 'gap'
    df = df.append(df2, ignore_index=True)
    df = df.sort_values(by=[0])
    df = df.rename(columns={0: 'to.id'})
    from_ids = df['to.id']
    from_ids = [0] + from_ids.tolist()[:-1]
    df.insert(0, 'from.id', from_ids)
    df.insert(0, 'chr', chrom)
    df.insert(2, 'from.coord', df['from.id']*res)
    df.insert(4, 'to.coord', df['to.id']*res)
    df.insert(6, 'size', df['to.coord']-df['from.coord'])
    if filter_func == signal_func.statFilter:
        for index, row in df.iterrows():
            if row['type'] == 'domain':
                df.loc[index, 'type'] = signal_func.add_boundary(p_values,
                                                                row['from.id'],
                                                                row['to.id'])
    bin_sig = pd.DataFrame(bin_sig)
    bin_sig.insert(0, 'chr', chrom)
    bin_sig.insert(1, 'from.coord', np.arange(len(bin_sig))*res)
    bin_sig.insert(2, 'to.coord', np.arange(1, len(bin_sig)+1)*res)
    df = df.rename(columns={0: 'mean.cf'})
    if bin_signal:
        return df, bin_sig
    return df


def find_not_gap(sig, w):
    n_bins = len(sig)
    gap = np.zeros(n_bins)
    for i in range(n_bins):
        if (sum(sig[i, max(1, i-w): min(i+w, n_bins)]) == 0):
            gap[i] = -0.5
    no_gaps = np.argwhere(gap != -0.5)
    no_gaps_list = []
    no_gaps = no_gaps.flatten()
    no_gaps += 1
    s = no_gaps[0]
    for i in range(1, len(no_gaps)):
        if no_gaps[i] != no_gaps[i-1]+1:
            e = no_gaps[i-1]
            no_gaps_list.append([s, e])
            s = no_gaps[i]
    no_gaps_list.append([s, no_gaps[-1]])
    return no_gaps_list


def binSignal(matrix, window, s):
    bin_sig = np.zeros(s)
    for i in range(1, s):
        lowerbound = max(0, i-window)
        upperbound = min(i+window, s)
        bin_sig[i] = np.mean(matrix[lowerbound:i, i:upperbound])
    return bin_sig


def plot_signal(r, peaks, r_sm=None):
    plt.clf()
    plt.plot(r, 'bo', markersize=4, alpha=0.5)
    if not isinstance(r_sm, None):
        plt.plot(r_sm, alpha=0.5)
    plt.plot(peaks, r[peaks], 'ro', markersize=4)
    plt.show()


def plot_matrix(file, peaks):
    plt.clf()
    plt.imshow(np.genfromtxt(file, delimiter='\t'), cmap='binary')
    # peaks = 150 +peaks
    plt.scatter(peaks, peaks)
    plt.show()


parser = argparse.ArgumentParser(description = "Calculate TopDop algorihm for a given input file.")
parser.add_argument("-i", "--input", help="input file (a .tsv contact map)", required=True)
parser.add_argument("-r", "--resolution", help="resolution (either 100k or 25k)", required=False)
parser.add_argument('-d', '--outputdir', help="directory to write output files to", required=False)
parser.add_argument('-p', '--peakfindfunc', help="Function used to finding peaks.", required=False)
parser.add_argument('-s', '--smoothfunc', help="Smoothing function.", required=False)
parser.add_argument('-f', '--filterfunc', help="Function used to filter outputs.", required=False)
parser.add_argument('-o', '--outputfile', help="Name of the output.", required=False)

args = vars(parser.parse_args())
if args["input"]:
	file = args["input"].strip()
	chrom = file.split(".")[-2].split("_")[-1]
if not args['resolution']:
	res = 100000
elif args['resolution'] in ["100k", "25k"]:
	res = int(args["resolution"][:-1])*1000
else:
	print("Wrong resolution value. Choose from [25k, 100k] and try again.")
	exit()
if args['outputdir']:
	outputdir = args['outputdir']
else:
	outputdir =  ''

if args['peakfindfunc']:
	if args['peakfindfunc']=="None": peakfindfunc = None
	elif args['peakfindfunc']=="find_min": peakfindfunc = peakfindfunc = signal_func.find_min
	elif args['peakfindfunc']=="detect_local_extrema": peakfindfunc = signal_func.detect_local_extrema
	elif args['peakfindfunc']=="find_peaks": peakfindfunc = signal_func.find_peaks
	elif args['peakfindfunc']=="find_peaks_2": peakfindfunc = signal_func.find_peaks_2
else: peakfindfunc =  ''

if args['smoothfunc']:
	if args['smoothfunc'] == "False": smoothfunc=False
	elif args['smoothfunc'] == "savgol_filter": smoothfunc=signal_func.savgol_filter
	elif args['smoothfunc'] == "smooth": smoothfunc=signal_func.smooth
	elif args['smoothfunc'] == "qspline": smoothfunc=signal_func.qspline
else: smoothfunc =  ''

if args['filterfunc']:
	if args['filterfunc']=="False": filterfunc = False
	elif args['filterfunc']=="statFilter": filterfunc = signal_func.statFilter
	elif args['filterfunc']=="filter_peaks": filterfunc = signal_func.filter_peaks
else: filterfunc =  ''

if args['outputfile']:
	outputfile = args['outputfile']
else: outputfile =  ''



df = TopDom(file=file, window=5, chrom=chrom, res=res,
            peak_find_funk=peakfindfunc,
            filter_func=filterfunc, smooth_func = smoothfunc,
            bin_signal=False)
#df.to_csv(outputdir+"dataframe_%s%s.csv" %("_"+outputfile+"__", chrom), index=False)
#bin_signal.to_csv(outputdir+"bin_signal_%s.csv" %chrom, index=False)
df.to_csv(outputdir+"domains_%s%s.csv" %("_"+outputfile+"_", chrom), index=False, columns=["chr", "from.coord", "to.coord", "type"])
