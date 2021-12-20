#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
import signal_func
import pandas as pd


def TopDom(file, window, chrom, res, peak_find_funk=None,
           smoove_func=False, filter_func=False, bin_signal=False):
    matrix = np.genfromtxt(file, delimiter='\t')
    s = matrix.shape
    assert s[0] == s[1], 'Unknwon Type of matrix file'
    n_bins = s[0]
    bin_sig = binSignal(matrix, window, n_bins)
    no_gap_region = find_not_gap(matrix, window)
    all_peaks = []
    df = []
    for i in no_gap_region:
        print('prosses:', i[0], i[1])
        ana_reg = bin_sig[i[0]:i[1]]
        if smoove_func is not False:
            ana_reg = smoove_func(ana_reg)
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
    df['type'] = 'domein'
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
            if row['type'] == 'domein':
                df.loc[index, 'type'] = signal_func.add_boundry(p_values,
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


file = 'contact_map.tsv'

df, bin_signal = TopDom(file=file, window=5, chrom='chr21', res=100000,
                        peak_find_funk=signal_func.detect_local_extrema,
                        filter_func=signal_func.statFilter, bin_signal=True)

df, bin_signal = TopDom(file=file, window=5, chrom='chr21', res=100000,
                        peak_find_funk=signal_func.find_min,
                        filter_func=signal_func.statFilter, bin_signal=True)

print(df)
print(bin_signal)
