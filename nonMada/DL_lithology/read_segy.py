# -*- coding: utf-8 -*-
"""
# Copyright (C) 2018 Guoyin Zhang and Yangkang Chen
# Python version: python 3 
#
# Reference:
# 
# Zhang, G., Z. Wang, and Y. Chen, 2018, Deep Learning for Seismic Lithology Prediction, Geophysical Journal International, 215, 1368–1387.
# https://github.com/agile-geoscience/xlines/blob/master/notebooks/
01_Read_and_write_SEG-Y.ipynb
"""

import numpy as np
import matplotlib.pyplot as plt
from obspy.io.segy.segy import _read_segy
import pandas as pd

########################### function defination  ##############################
def segy_inf(segy):
    '''segy: data of SEGYFile format'''
    inline = np.array([t.header.trace_sequence_number_within_line 
                       for t in segy.traces])
    xline = np.array([t.header.original_field_record_number 
                      for t in segy.traces])
    coor_x = np.array([t.header.source_coordinate_x for t in segy.traces])
    coor_y = np.array([t.header.source_coordinate_y for t in segy.traces])
    print('inline range: %.f - %.f' %(min(inline),max(inline)))
    print('xline range: %.f - %.f' %(min(xline),max(xline)))
    print('Two-point coordinates: (%.f, %.f), (%.f, %.f)' % (min(coor_x), 
                                   min(coor_y),max(coor_x),max(coor_y)))
    '''plot the survey'''
    #plt.plot(xline, inline, '.')
    #plt.plot(coor_x, coor_y, '.')


def extract_line(segy, line_num):
    '''extract the exact line in the segy file
    segy: data of SEGYFile format
    line_num: number of line from segy'''
    data = np.stack(t.data for t in segy.traces 
                    if t.header.trace_sequence_number_within_line == line_num)
    data = data.T
    return data


def extract_line_v2(segy, inline, xline_min, xline_max, 
                    time_start, time_stop, time_step):
    '''extract line with time gate and trace gate'''
    line = pd.DataFrame()
    for xx in np.arange(xline_min, (xline_max+1)):  #2080
        if xx % 50 == 0: print("==== extracting xline %d" %xx)
        trace = extract_seismic(segy, line=inline, xline=xx, 
                                  time_start=time_start, time_stop=time_stop, 
                                  time_step=time_step, name='amp')
        line = pd.concat([line, trace.amp], axis=1)
    return line


def extract_trace(segy, line_num, xline_num):
    '''extract the exact trace in the segy file according to line and xline num
    segy: data of SEGYFile format
    line_num: number of line from segy
    xline_num: number of the xline from segy'''
    inline = np.array([t.header.trace_sequence_number_within_line 
                       for t in segy.traces])
    xline = np.array([t.header.original_field_record_number 
                      for t in segy.traces])
    num = (line_num-min(inline)) * (max(xline)-min(xline)+1) + \
          (xline_num-min(xline))
    trace = segy.traces[num].data
    #another way to extract the data directly
    #trace = np.stack(t.data for t in segy.traces 
    #                if t.header.trace_sequence_number_within_line == line_num 
    #                and t.header.original_field_record_number == xline_num)
    return trace
    

def extract_seismic(segy, line, xline, 
                    time_start, time_stop, time_step, 
                    name='seismic'):
    '''extract the seismic data near well in a time period'''
    '''output data type: DataFrame
    
    segy: imported segy data
    line: line number to be extracted
    xline: xline number to be extracted
    time_start: start time to be extracted
    time_stop: stop time to be extracted
    time_step: seismic_interval in segy file
    name: name of the data, like segy, ampulitute, rms
    '''
    trace = extract_trace(segy, line, xline)
    original_time = segy.traces[0].header.delay_recording_time
    
    index_start = int((time_start-original_time)/time_step)
    index_stop = int((time_stop-original_time)/time_step)
    trace_target =trace[index_start:index_stop]
    
    trace_target = pd.DataFrame(np.array(trace_target), columns=[name])
    time = np.arange(time_start, time_stop, time_step)
    trace_target = pd.concat([pd.DataFrame(time, columns=['time']), 
                              trace_target], axis=1)
    return trace_target

        
def plot_line(data, xleft, xright, yup, ydown, cmap="RdBu", percent=99):
    vm = np.percentile(data, percent)
    plt.figure(figsize=(16,8))
    plt.imshow(data, cmap=cmap, interpolation='bilinear', 
               vmin=-vm, vmax=vm, 
               extent=(xleft, xright, ydown, yup),
               aspect='auto')  #'nearest' 'bilinear' 
    plt.colorbar()
    plt.tight_layout()
    #plt.xlim(xmin=400, xmax=800)
    #plt.ylim(ymax=800,ymin=1200)
    #plt.savefig("lin688_original_seis2.png", dpi=250)
    

def plot_line_attribute(data, cmap="RdBu", percent=99):
    vm = np.percentile(data, percent)
    plt.figure(figsize=(16,8))
    plt.imshow(data, cmap=cmap, interpolation='bilinear', 
               vmin=0, vmax=vm, aspect='auto')
    plt.colorbar(label="Amplitude", shrink=0.5)
    #plt.xlim(xmin=400, xmax=800)
    #plt.ylim(ymax=800,ymin=1200)
   
#%%
########################### main process  ##############################
if __name__ == "__main__":
    
    '''return a SEGYFile object '''
    segy = _read_segy('Data/CX_pst_3710_3720.sgy', headonly=True)
    
    '''read the header and trace'''
    #print(segy.binary_file_header)
    #print(segy.traces[492].header)
    segy_inf(segy)
    
    #print(segy.textual_file_header.decode(encoding='cp037'))
    x = np.array(list(segy.textual_file_header.decode(encoding='cp037')))
    print('\n'.join(''.join(row) for row in x.reshape((40, 80)))) 
    print(segy.binary_file_header)
    print(segy.traces[2].header)

    '''etract line and plot the line'''
    line3714 = extract_line(segy, 3714)   #rows are traces
    plot_line_attribute(line3714, cmap="seismic", percent=99)
    
    #%%
    '''calculate seismic attribute'''
    import bruges
    segy = _read_segy('Data/CX_pst_3710_3720.sgy', headonly=True)
    data = np.stack(t.data for t in segy.traces)
    print(segy.textual_file_header.decode())
    print(segy.binary_file_header)
    
    dt = segy.traces[0].header.sample_interval_in_ms_for_this_trace / 1e6
    
    #energy
    energy = bruges.attribute.energy(segy, duration=0.025, dt=dt)
    plot_line_attribute(energy, cmap="jet")  #cmap: "viridis" "Greys"
    #plt.savefig("lin688_energy.png", dpi=250)
    
    #similarity
    similarity = bruges.attribute.similarity(line3714, duration=0.025, dt=dt, 
                                             step_out=1, lag=0)
    plot_line_attribute(similarity, cmap="viridis_r", percent=99.9)
    #plt.savefig("lin688_similarity2.png", dpi=250)
    
    #trace = extract_trace(segy, 688, 390).data   



    #%%
    ################################### Write data  ###########################
    #write
    from obspy.core import AttribDict
    from obspy.core import Stats
    
    from obspy.core import Trace, Stream
    from obspy.io.segy.segy import SEGYBinaryFileHeader
    from obspy.io.segy.segy import SEGYTraceHeader
    
    out = Stream(Trace(t, header=dict(delta=dt)) for t in data)
    out.stats = Stats(dict(textual_file_header=segy.textual_file_header))
#    out.stats.textual_file_header += """energy volume.
#    Generated: 18 Sep 2016 by Matt Hall matt@agilegeoscience.com.
#    Algorithm: github.com/agile-geoscience/bruges.attribute.similarity.
#    Parameters: duration=0.16 s, dt=0.002 s.""".encode('utf-8')
    
    out.write('out1.sgy', format='SEGY',data_encoding=1)
                                        #encode 1 for IBM,5 for IEEE
    #%%
    #####################  Validate the output data  ##########################
    s = _read_segy('out1.sgy', headonly=True)
    s.textual_file_header[:15*80]  # First 15 lines.
    s.traces[0].header
    segy.traces[0].data.shape
    d=np.arange(751)
    
    segy = _read_segy('Data/3X_75_PR.sgy', headonly=True)
    segy.write('file2.segy') 
    segyf = _read_segy('file2.segy', headonly=True)
    segyf.traces[10].header
    
    plt.figure(figsize=(16,8))
    plt.imshow(np.vstack(t.data for t in s.traces).T,
               cmap='viridis',
               aspect='auto')
    plt.colorbar()
    plt.show()


























