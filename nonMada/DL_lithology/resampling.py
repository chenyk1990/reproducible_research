# -*- coding: utf-8 -*-
"""
# Copyright (C) 2018 Guoyin Zhang and Yangkang Chen
# Python version: python 3 
#
# Reference:
# 
# Zhang, G., Z. Wang, and Y. Chen, 2018, Deep Learning for Seismic Lithology Prediction, Geophysical Journal International, 215, 1368–1387.
"""

from scipy import interpolate  
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

#####################  function definition    #################################

def depth_to_time(vsp_depth, vsp_time, depth, kind="quadratic", plot=False):
    '''
    This function is use to convert well depth to time, according to 
    depth-time relationship like vsp. 
    
    The range of vsp_depth must be larger than depth
    
    vsp_depth: depth in the original time-depth relationship
    vsp_time: time in the original time-depth relationship
    depth: depth for converting, DataFrame
    kind: "nearest","zero","slinear","quadratic","cubic"
    plot: if the default plot is True, the figure will be plotted
    '''    
    f = interpolate.interp1d(vsp_depth, vsp_time, kind=kind)
    time = f(depth)
    depth_time = pd.concat([depth, 
                            pd.DataFrame(time, columns=['time'])], 
                            axis=1)    
#    depth_time = pd.concat([pd.DataFrame(depth, columns=['depth']), 
#                            pd.DataFrame(time, columns=['time'])], 
#                            axis=1)
    if plot==True:
        plt.plot(vsp_depth, vsp_time, "b--", label='original')
        plt.plot(depth, time, "r-", linewidth = 3, label='converted') 
        plt.legend(loc="lower right") 
        plt.xlabel(u"Depth(m)")
        plt.ylabel(u"Time(ms)")
        plt.title(u"Time conversion")        
        plt.show()
    return depth_time


def resample(x, y, time_start_re, time_stop_re, time_step_re, 
             kind="quadratic", plot=False, logname='log'):
    '''
    This function is use to resample the tagert variable to a new space.
    output data type: pd.DataFrame
    
    x: depth in the original time-depth relationship
    y: time in the original time-depth relationship
    kind: "nearest" and "zero" are better for discrete variables
          "slinear","quadratic", and "cubic" are better for continous variables 
    plot: if the default plot is True, the figure will be plotted
      
    time_start_re:  time_start_resampling
    time_stop_re:
    time_step_re:    # ms
    '''
    time_resampling = np.arange(time_start_re, time_stop_re, time_step_re)
    
    f = interpolate.interp1d(x, y, kind=kind)
    y_r = f(time_resampling)
    xy_r = pd.concat([pd.DataFrame(time_resampling, columns=['time']), 
                            pd.DataFrame(y_r, columns=[logname])], 
                            axis=1)
    if plot==True:
        plt.figure(figsize=(12,3))
        plt.plot(x, y, "b--", label='original')
        plt.plot(time_resampling, y_r, "ro", label='resampled') 
        plt.xlabel(u"Time(ms)")
        plt.ylabel(u"Value")
        plt.title(u"Data resampling")
        plt.legend(loc="lower right")  
        plt.show()
    return xy_r   
