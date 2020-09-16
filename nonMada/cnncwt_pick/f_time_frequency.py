# -*- coding: utf-8 -*-
"""
# Copyright (C) 2019 Guoyin Zhang and Yangkang Chen
# Python version: python 3 
#
# Reference:
# Chen, Y., G. Zhang, M. Bai, S. Zu, Z. Guan, and M. Zhang, 2019, Automatic Wave- form Classification and Arrival Picking Based on Convolutional Neural Network, Earth and Space Science, 6, 1244-1261.
# Zhang, G., C. Lin, and Y. Chen, 2020, Convolutional Neural Networks for Microseismic Waveform Classification and Arrival Picking, Geophysics, 85, WA227–WA240.
# Zhang, G., Z. Wang, and Y. Chen, 2018, Deep Learning for Seismic Lithology Prediction, Geophysical Journal International, 215, 1368–1387.
"""

import matplotlib.pyplot as plt
import numpy as np
import pywt 
from obspy.io.segy.segy import _read_segy
import scipy.signal as signal

    
def spectrum_fft_filter(time, data, fft_size, 
                        freq_set=False, freq_min=0, freq_max=200, 
                        withfilter=False, filter_slop=10, filter_value=100, 
                        btype='low', analog=False, output='ba'):
    '''plot the spectrum of one-dimensional data, filter can be applied
    
    time: time of data, ms
    data: the targeted one-dimensional data
    fft_size: the window size for fft 
    
    freq_set: False or True
              if True, freq_min and freq_max can be set for x axes
    freq_min: minimum of frequency for x axes
    freq_max: maximum of frequency for x axes
    
    withfilter: False or True
                if true, filter will be applied to data
    filter_slop: the slope of the filter
    filter_value: the maximum or minimum value that you want to be filtered
    btype: {‘lowpass’, ‘highpass’, ‘bandpass’, ‘bandstop’}, optional
           The type of filter. Default is ‘lowpass’.
    analog : bool, optional
             When True, return an analog filter, otherwise a digital filter is\
             returned.
    output : {‘ba’, ‘zpk’, ‘sos’}, optional
             Type of output: numerator/denominator (‘ba’), pole-zero (‘zpk’), \
             or second-order sections (‘sos’). Default is ‘ba’.
    '''
    x = np.array(time).reshape(len(time),)
    y = np.array(data).reshape(len(data),)
    
    sampling_rate = 1000 * (len(x)-1)/(max(x)-min(x))
    
    '''one-dimensional discrete Fourier Transform for real input'''
    magnitude = np.fft.rfft(y, fft_size) # 幅度
    amplitude = np.abs(magnitude)*2/fft_size #绝对幅值或振幅
    amplitude[0] = 0 #直流分量设为0
    
    '''filter'''
    if withfilter == True:
        b, a = signal.butter(filter_slop, filter_value/(sampling_rate/2), 
                             btype=btype, analog=analog, output=output)
        y_f = signal.filtfilt(b, a, y)
        magnitude_f = np.fft.rfft(y_f, fft_size) # 幅度
        amplitude_f = np.abs(magnitude_f)*2/fft_size #绝对幅值或振幅
        amplitude_f[0] = 0 #直流分量设为0
    
    '''frequency distribution'''
    df = sampling_rate/ fft_size
    freqs = np.linspace(0, fft_size/2, fft_size/2+1) * df
    
    '''plot the map'''
    plt.figure(figsize=(10,6))
    plt.subplot(211)
    plt.plot(x, y, 'k') #'steelblue'
    plt.fill_between(x, y, y2=0, where=y >= 0, facecolor='k', interpolate=True)
    if withfilter == True:
        plt.plot(x, y_f, 'r')
    plt.xlabel(u"Time(ms)")
    plt.ylabel(u"Value")
    plt.ylim(-5,5)
    plt.xlim(940, 1375)
    plt.title(u"Data and Spectrum")
    
    plt.subplot(212)
    plt.vlines(freqs, [0], amplitude, 'r') #'steelblue'  'darkred'
    if withfilter == True:
        plt.vlines(freqs, [0], amplitude_f, 'r')
    plt.xlabel(u"Frequency(Hz)")
    plt.ylabel(u"Amplitude")
    if freq_set == True:
        plt.xlim(xmax=freq_max, xmin=freq_min)
    plt.subplots_adjust(hspace=0.4)
    plt.tight_layout()
    plt.show() 
    if withfilter == True:
        return y_f
    return freqs,amplitude


def spectrum_cwt(t, data, wavelet='mexh', widths=200, cmap='RdBu', 
                 colorscale=1, contour_dvision=41, freqmax=100,
                 figsize=(12,3), plot=True,ymin=5, ymax=30, vmax=None):
    '''
    cwt: contineous wavelet transforms for seismic
    '''
    '''
    t: in ms
    data: 1d array
    wavelet: 
        ['haar', 'db', 'sym', 'coif', 'bior', 'rbio', 'dmey', 'gaus', 'mexh', 
        'morl', 'cgau', 'shan', 'fbsp', 'cmor'] short for:
        ['Haar', 'Daubechies', 'Symlets', 'Coiflets', 'Biorthogonal', 
        'Reverse biorthogonal','Discrete Meyer(FIR Approximation)','Gaussian', 
        'Mexican hat wavelet', 'Morlet wavelet', 'Complex Gaussian wavelets', 
        'Shannon wavelets', 'Frequency B-Spline wavelets', 
        'Complex Morlet wavelets']
    widths: the number of frequencies
    colorscale: change color scale, the smaller the darker the color
    contour_division: the number of color divided
    freqmax: the maximum of frequency
    plot: True or False
    '''
    t = np.array(t).reshape(len(t),)
    data = np.array(data).reshape(len(data),)
    
    widths = np.arange(1, widths)
    sampling_rate = 1000 * len(t)/(max(t) - min(t))
    cwtmatr, freqs = pywt.cwt(data, widths, wavelet, 1/sampling_rate)
    
    maximum = colorscale * max(abs(cwtmatr.max()), abs(cwtmatr.min()))
    if vmax != None:  maximum = vmax
    
    cwtmatr_c = np.clip(cwtmatr, -maximum, maximum, cwtmatr)
    levels = np.linspace(-maximum, maximum, contour_dvision)
    
    if plot==True:
        plt.figure(figsize=figsize)
        plt.contourf(t, freqs, cwtmatr_c, levels=levels, cmap=cmap)
        #plt.ylim(min(freqs), freqmax) #
        #plt.yscale('log')
        plt.ylim(ymin, ymax)
        #plt.xlim(t.min(), t.max()+1)
        plt.xlabel(u"Time(ms)")
        plt.ylabel(u"Frequency(Hz)")
        plt.title(u"Time-frequency Spectrum")
        plt.colorbar()
        
        plt.tight_layout() 
        plt.show()
    
    
    return cwtmatr, freqs


def spectrum_cwt2(t, data, wavelet='mexh', widths=200, cmap='RdBu', 
                 colorscale=1, contour_dvision=41, freqmax=100,
                 figsize=(4.5,4), plot=True): #(3.5,8)
    '''
    similiar to spectrum_cwt
    exchange the position of X and Y of cwt map
    '''
    t = np.array(t).reshape(len(t),)
    data = np.array(data).reshape(len(data),)
    
    widths = np.arange(1, widths)
    sampling_rate = 1000*len(t)/(max(t) - min(t))
    cwtmatr, freqs = pywt.cwt(data, widths, wavelet, 1/sampling_rate)
     
    #maximum = colorscale*max(abs(cwtmatr.max()), abs(cwtmatr.min()))
    maximum = colorscale*4
    cwtmatr_c = np.clip(cwtmatr, -maximum, maximum, cwtmatr)
    levels = np.linspace(-maximum, maximum, contour_dvision)
    if plot==True:
        plt.figure(figsize=figsize)
        plt.contourf(freqs, t, cwtmatr_c.T, levels=levels, cmap=cmap)
        #plt.ylim(min(freqs), freqmax) #
        plt.xlim(10, 80)
        plt.ylim(1055, 1096) #(1055, 1096) (940, 1375)  
        plt.xlabel(u"Frequency(Hz)")
        plt.ylabel(u"Time(ms)")
        #plt.title(u"Time-frequency Spectrum")
        plt.colorbar()
        ax = plt.gca()
        ax.invert_yaxis()
        plt.tight_layout() 
        plt.show()
    return cwtmatr, freqs

"""
#####################       Main processes    #################################
if __name__ == "__main__":
    
    import read_segy as rsegy
    import resampling as rs

    segy = _read_segy('data/seis/CX_ZJ_ori.SGY', headonly=True)
#    time_start = 900.0
#    time_stop = 1500.0
#    time_step = 2.0 # ms
#    inline = 3662
#    xline = 1938
    trace = rsegy.extract_seismic(segy, line=3662, xline=1938, 
                time_start = 810, time_stop= 1500, time_step=2, name='amp')
    trace_res = rs.resample(trace.time, trace.amp, 
                                    time_start_re = 820,
                                    time_stop_re = 1400,
                                    time_step_re = 1,
                                    kind="quadratic", plot=False, 
                                    logname='amp')
    cc,ff = spectrum_cwt(trace_res.time, trace_res.amp/10000, 
                    colorscale=1, widths=100, wavelet='morl', freqmax=100)
    freqs,amplitude = spectrum_fft_filter(trace_res.time, trace_res.amp/10000, 
                        fft_size = 2000, 
                    freq_set=True, freq_min=0, freq_max=100, 
                    withfilter=False, filter_value=60, filter_slop=10,  
                    btype='low', analog=False, output='ba')

    plt.figure(figsize=(3,8))
    plt.plot(trace_res.amp/10000, trace_res.time, 'k') #'steelblue'

    plt.xlabel(u"Amplitude")
    plt.ylabel(u"Time(ms)")
    plt.ylim(940, 1375)
    plt.xlim(-5,5)
    plt.title(u"Seismic")
    plt.fill_betweenx(trace_res.time, trace_res.amp/10000, x2=0, 
                      where=trace_res.amp/10000 >= 0, 
                      facecolor='k', interpolate=True)
    ax = plt.gca()
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.show()

"""








