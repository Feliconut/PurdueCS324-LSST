# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 19:58:53 2021

@author: Ankita Mishra
"""
# Refactored 3/23/2021 by Xiaoyu Liu

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


#Find peaks
def find_maxima(x, y):
    peaks = find_peaks(y, height=1, threshold=1, distance=1)
    peak_height = peaks[1]['peak_heights']  #list of the heights of the peaks
    peak_pos = x[peaks[0]]  #list of the peaks positions
    print(peak_pos)
    print(peak_height)
    return peak_pos, peak_height


#Finding the minima
def find_minima(x, y):
    y2 = -y
    minima = find_peaks(y2)
    min_pos = x[minima[0]]  #list of the minima positions
    min_height = y2[minima[0]]  #list of the mirrored minima heights
    return min_pos, -min_height


def visualize_max_min(x, y):
    #Plotting
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(x, y)
    ax.scatter(*find_maxima(x, y), color='r', s=15, marker='D', label='Maxima')
    ax.scatter(*find_minima(x, y),
               color='gold',
               s=15,
               marker='X',
               label='Minima')
    ax.legend()
    ax.grid()
    plt.show()


if __name__ == '__main__':
    #defining the x and y arrays
    x = np.linspace(0, 10, 100)
    y = x * np.random.randn(100)**2
    visualize_max_min(x, y)
