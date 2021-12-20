# -*- coding: utf-8 -*-

import scipy.signal as scpsi
import numpy as np
from sklearn.preprocessing import scale
from scipy.stats import mannwhitneyu

'''smoothing'''


def savgol_filter(r):
    # slabe
    r_sm = scpsi.savgol_filter(r, 21, 5)
    return r_sm


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def qspline(r):
    time = np.arange(len(r))
    filtered = scpsi.qspline1d_eval(scpsi.qspline1d(r), time)
    return filtered


'''finding min'''


def find_min(r):
    peaks = scpsi.argrelextrema(r, np.less)[0]
    return peaks


def find_peaks(r, d=None):
    peaks, _ = scpsi.find_peaks(1-r, height=0, distance=d)
    return peaks


def find_peaks_2(r):
    peaks = scpsi.find_peaks_cwt(1-r, np.arange(1, 10))
    return peaks


'''filtering'''


def filter_peaks(r, peaks, t=0.02):
    prominences = scpsi.peak_prominences(1-r, peaks)[0]
    p_t = prominences > t
    peaks = peaks[p_t]
    return peaks


'''original functions'''


def detect_local_extrema(x):
    x = np.append(x, [0])
    local_min = []
    local_max = []
    if len(x) <= 3:
        local_min = np.argmin(x)
        local_max = np.argmax(x)
        return local_min, local_max
    else:
        x = data_norm(x)
        cp = change_point(x)
        if len(cp) <= 2 or len(cp) == len(x):
            return local_min  # local_max
        for i in range(1, (len(cp)-1)):
            if x[cp[i]] >= x[cp[i]-1] and x[cp[i]] >= x[cp[i]+1] :
                local_max.append(cp[i])
            if x[cp[i]] < x[cp[i]-1] and x[cp[i]] < x[cp[i]+1]:
                local_min.append(cp[i])
            min_val = min(x[cp[i-1]], x[cp[i]])
            max_val = max(x[ cp[i-1] ], x[ cp[i] ] )
            if(min(x[cp[i-1]:cp[i]]) < min_val):
                a = cp[i-1] - 1 + np.argmin(x[cp[i-1]:cp[i]])+1
                local_min.append(a)
            if max(x[cp[i-1]:cp[i]] ) > max_val:
                b=  cp[i-1] - 1 + np.argmax(x[cp[i-1]:cp[i]] )
                local_max.append(b)
    return np.sort(local_min)  # local_max


def data_norm(y):
    y_diff = np.diff(y)
    scale_y = 1/np.mean(np.abs(y_diff))
    for i in range(1, len(y)):
        y[i] = y[i-1] + (y_diff[i-1]*scale_y)
    return y


def change_point(y):
    x = np.arange(1, len(y)+1)
    cp, i = [0], 0
    Fv = [None] * len(y)
    # Ev = [None] * len(y)
    Fv[0] = 0
    while i < len(y)-1:
        j = i+1
        Fv[j] = np.sqrt((x[j]-x[i])**2 + (y[j]-y[i])**2)
        while j < len(y)-1:
            j += 1
            # Ev[j] = np.sum(np.abs((y[j]-y[i])*x[i+1:j-1] - (x[j]-x[i])*y[i+1:j-1] -(x[i]*y[j]) + (x[j]*y[i]))) / np.sqrt((x[j]-x[i])**2 + (y[j]-y[i])**2)
            P = np.sqrt((x[j]-x[i])**2 + (y[j] - y[i])**2) - (np.sum(np.abs((y[j]-y[i])*x[i+1:j] - (x[j]-x[i])*y[i+1:j]- (x[i]*y[j]) + (x[j]*y[i])))/ np.sqrt((x[j]-x[i])**2 + (y[j] - y[i])**2))
            Fv[j] = P
            if(Fv[j] == None or Fv[j] == None):
                j = j-1
                cp = cp + [j]
                break
            if Fv[j] < Fv[j-1]:
                j = j-1
                cp = cp + [j]
                break
        i = j
    cp = cp + [len(y)-1]
    return cp


def statFilter(matrix, w, no_gap_region, peaks):
    n = matrix.shape[0]
    p_values = np.ones(n)
    scale_matrix = matrix
    for i in range(1,2*w):
        scale_matrix[np.arange(0, n-i), np.arange(i, n)] = scale(matrix.diagonal(i))
    for i in no_gap_region:
        p_values[i[0]-1:i[1]] = get_pvalues(scale_matrix[i[0]-1:i[1], i[0]-1:i[1]], w)
    good_pos = np.where(p_values < 0.05)[0]
    good_pos +=1
    good_peaks = np.intersect1d(good_pos, peaks)
    return good_peaks, p_values



def get_pvalues(matrix, size, scale=1):
    n_bins = matrix.shape[0]
    p_vals = np.ones(n_bins)
    for ind, i in enumerate(range(1, n_bins)):
        dia = get_diamond_matrix(matrix, i, size)
        dia = dia.flatten()
        ups = get_upstream_triangle(matrix, i, size=size)
        downs = get_downstream_triangle(matrix, i, size=size)
        p = np.concatenate((ups, downs), axis=None)
        w, p_val = mannwhitneyu(dia, p, alternative='less')
        p_vals[ind] = p_val
    return p_vals

def get_diamond_matrix(matrix, i, size):
    lower = max(0, i-size)
    upper = min(matrix.shape[1], i+size)
    return matrix[lower:i, i:upper]

def get_upstream_triangle(matrix, i, size):
    lower = max(0, i-size-1)
    s_matrix = matrix[lower:i, lower:i]
    s = s_matrix.shape[0]
    idu = np.triu_indices(s, k=1)
    return s_matrix[idu]

def get_downstream_triangle(matrix, i, size):
    upper  = min(i+size, matrix.shape[1])
    s_matrix = matrix[i:upper, i:upper]
    s = s_matrix.shape[0]
    idu = np.triu_indices(s, k=1)
    return s_matrix[idu]


def add_boundry(p_vals, s, e):
    good_pval = np.where(p_vals[s:e] < 0.05, 0, 1)
    if np.sum(good_pval) == 0:
        return 'boundary'
    else:
        return 'domein'
    
    
    
# # #%% savgol jest slaby
# r_sm = savgol_filter(r)
# peaks1 = find_min(r_sm)
# plot_signal(r, peaks1, r_sm=r_sm)
# plot_matrix(file, peaks1)

# #%% sparwdza sie trzeba ustawic d
# peaks2 = find_peaks(r)
# plot_signal(r, peaks2)
# plot_matrix(file, peaks2)

# #%% nie jest zle
# peaks3 = find_peaks_2(r)
# plot_signal(r, peaks3)
# plot_matrix(file, peaks3)

# #%% trzeba by porownac z org
# r_sm = smooth(r, 3)
# peaks4 = find_min(r_sm)
# plot_signal(r, peaks4, r_sm=r_sm)
# plot_matrix(file, peaks4)

# #%% smooth pomaga
# peaks5 = find_min(r)
# plot_signal(r, peaks5)
# plot_matrix(file, peaks5)

# #%% ograniczanie wielkosci pikow przez d
# r_sm = smooth(r, 3)
# peaks6 = find_peaks(r_sm, d =10)
# plot_signal(r, peaks6, r_sm=r_sm)
# plot_matrix(file, peaks6)


# #%% nie wiele wnosi

# r_sm = qspline(r)
# peaks7 = find_peaks(r_sm)
# plot_signal(r, peaks7, r_sm=r_sm)
# plot_matrix(file, peaks7)


# #%%
# peaks8 = find_peaks(r)
# peaks8 = filter_peaks(r, peaks8, 0.03)
# plot_signal(r, peaks8, r_sm=r_sm)
# plot_matrix(file, peaks8)
