# -*- coding: utf-8 -*-
"""
# Copyright (C) 2018 Guoyin Zhang and Yangkang Chen
# Python version: python 3 
#
# Reference:
# 
# Zhang, G., Z. Wang, and Y. Chen, 2018, Deep Learning for Seismic Lithology Prediction, Geophysical Journal International, 215, 1368–1387.
"""

#%% ##################      Pakages imported   ################################
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib import gridspec

import numpy as np
import pandas as pd
from obspy.io.segy.segy import _read_segy
from keras.optimizers import Adam, SGD
from keras.models import load_model
from keras.utils import np_utils
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import f1_score


from keras.models import Sequential, Model
from keras.layers import Dense, Activation,normalization, Dropout
from keras.layers import Conv1D, MaxPooling1D, Conv2D, MaxPooling2D, Flatten
from keras.optimizers import SGD,RMSprop,Adam
from keras import backend as K
from keras import callbacks

import read_segy as rsegy
import read_las as rlas
import resampling as rs
import time_frequency as tf

#%% ##################      Function  defination   ############################
class data:
    '''
    generating samples for one trace or one well
    '''
    def cwt2d(sample, timename, yname, xname, num=10, interval=1,
              figsize=(5,5), plot=False):
        '''
        Generating cwt map samples for CWT-CNNs
        '''
        '''
        sample: integrated logs and seismic data file, DataFrame
        timename: the name of the time column in sample, string
        yname: the name of the target parameter, string
        xname: the name of the input parameter, string
        num: number of data beside the targeted data, 
                    total num will be 2*sample+1
        interval: time interval between samples
        figseze: the figure size of the cwtmap
        plot: True or False, plot cwt map or not
        ''' 
        length = sample.shape[0]
        time = sample.loc[num*interval: (length-num*interval-1), timename]
        y = sample.loc[num*interval: (length-num*interval-1), yname]
        seis = sample.loc[num*interval: (length-num*interval-1), xname]    
        
        data = np.zeros((int(length-num*interval*2),
                         len(xname), 
                         15, 
                         int(2*num+1)))
                
        for m in np.arange(0,len(xname)):
            '''cwt operation'''
            cwtmatr, freqs = tf.spectrum_cwt(
                            sample.loc[:, timename], 
                            sample.loc[:, xname[m]], 
                            colorscale=1, wavelet='morl', widths=100, 
                            freqmax=100, figsize=figsize, plot=plot)
            '''cwt map sampling'''
            for i in np.arange(num*interval, length-num*interval):
                data[(i-num*interval),m,:,:] = cwtmatr[
                 [9,10,11,12,13,14,15,17,19,22,26,31,39,53,80],:][:,
                 list(np.arange((i-num*interval),(i+num*interval+1),interval))]
        y = np.array(y)
        time = pd.DataFrame(time)
        '''
        data: cwt maps
        y: the target parameter
        time: the time of the samples
        seis: the original sesimic values
        freqs: the resampled frequencies
        '''
        return data, y, time, seis, freqs


    def cwt2d_seis(sample, timename, xname, num=10, interval=1,
                  figsize=(5,5), plot=False):
        '''
        Generating cwt map samples for CWT-CNN line processing.
        Similiar with cwt2d, but without y in the result
        '''
        length = sample.shape[0]        
        data = np.zeros((int(length-num*interval*2),
                         len(xname), 15, int(2*num+1)))
        for m in np.arange(0,len(xname)):
            cwtmatr, freqs = tf.spectrum_cwt(
                            sample.loc[:, timename], 
                            sample.loc[:, xname[m]], 
                            colorscale=1, wavelet='morl', widths=100, 
                            freqmax=100, figsize=figsize, plot=plot)
            for i in np.arange(num*interval, length-num*interval):
                data[(i-num*interval),m,:,:] = cwtmatr[
                 [9,10,11,12,13,14,15,17,19,22,26,31,39,53,80],:][:,
                 list(np.arange((i-num*interval),(i+num*interval+1),interval))]
        return data

    
    def raw1d(sample, timename, yname, xname, num=10, interval=1):
        '''
        Generating 1d vector samples for DNNs, CNNs
        '''
        '''        
        sample: integrated logs and seismic data file, DataFrame
        timename: the name of the time column in sample, string
        yname: the name of the target parameter, string
        xname: the name of the input parameter, string
        num: number of data beside the targeted data, 
                    total num will be 2*sample+1
        interval: time interval between samples
        '''
        length = sample.shape[0] 
        data = pd.DataFrame()
        for i in np.arange(num*interval, length-num*interval):
            #选择第i个目标
            data_t = pd.DataFrame()
            data_t = pd.concat([data_t, sample.loc[i, [timename]]])
            data_t = pd.concat([data_t, sample.loc[i, [yname]]])
            for k in xname:
                #选择变量
                for j in np.arange(-num*interval, num*interval+1, interval):
                    #根据取样间隔和长度选择特定序号变量
                    data_t = pd.concat([data_t, sample.loc[i+j,[k]]])
            data_t = pd.DataFrame(np.array(data_t.T))
            data = pd.concat([data, data_t])
        
        '''generate the names for the columns'''
        name = [timename]
        name.append(yname)
        for k in xname:
            for j in np.arange(1, num*2+2):
                name_j ="%s%.0f" %(k, j)
                name.append(name_j)
        '''rename the index and columns of the data'''
        data = pd.DataFrame(np.array(data), columns=name)
        return data


    def raw1d_seis(sample, timename, xname, num=10, interval=1):
        '''
        Generating 1d vector samples for DNN, CNN line processing.
        Similar with raw1d.
        '''
        length = sample.shape[0] 
        data = pd.DataFrame()
        for i in np.arange(num*interval, length-num*interval):
            #选择第i个目标
            data_t = pd.DataFrame()
            for k in xname:
                #选择变量
                for j in np.arange(-num*interval, num*interval+1, interval):
                    #根据取样间隔和长度选择特定序号变量
                    data_t = pd.concat([data_t, sample.loc[i+j,[k]]])
            data_t = pd.DataFrame(np.array(data_t.T))
            data = pd.concat([data, data_t])
        '''rename the index and columns of the data'''
        data = pd.DataFrame(np.array(data))
        return data


class dataset:
    '''
    generating data set for multiple traces or wells
    '''
    def cwt2d(wells, folderpath='data/zj_las', num=10, interval=2, 
              figsize=(8,5), cwtplot=True):
        '''
        Generating cwt map data set for CWT-CNNs
        '''
        '''
        wells: a list of well names for sampling
        folderpath: the folder path for log files
        
        num: number of data beside the targeted data, 
                    total num will be 2*sample+1
        interval: time interval between samples
        figsize: the figure size of the cwtmap
        cwtplot: True or False, plot cwt map or not
        '''
        for i in np.arange(0, len(wells)):
            ## read the las log file
            print("==== Generating the dataset of well %s"%wells[i])
            print("==== ==== resampling")
            logs = rlas.las_to_df(file='%s/%s' %(folderpath,wells[i]), 
                                  interval=False,
                                  depth_start=0, depth_stop=0, nan=-999.25)
            
            ## read the time gate and trace location from inf file
            inline = inf.loc[inf.name==wells[i], ['inline']].iloc[0,0]
            xline = inf.loc[inf.name==wells[i], ['xline']].iloc[0,0]
            time_start = inf.loc[inf.name==wells[i], ['start']].iloc[0,0]
            time_stop = inf.loc[inf.name==wells[i], ['stop']].iloc[0,0]
            time_step = 1  #resampling time steps
            
            ## read and resample the log data
            lith_res = rs.resample(logs.TIME, logs.LITH, 
                                   time_start_re = time_start+c,
                                   time_stop_re = time_stop-c,
                                   time_step_re = time_step,
                                   kind="nearest", plot=False, logname='lith')
            ## read and resample the seismic data
            trace = rsegy.extract_seismic(segy, line=inline, xline=xline, 
                                      time_start=800, time_stop=1500, 
                                      time_step=2, name='amp')
            trace_res = rs.resample(trace.time, trace.amp, 
                                    time_start_re = time_start+c,
                                    time_stop_re = time_stop-c,
                                    time_step_re = time_step,
                                    kind="quadratic", plot=False, 
                                    logname='amp')
            sample = pd.concat([lith_res, trace_res.amp], axis=1)
            
#            '''read and resample the rms data'''
#            trace = rsegy.extract_seismic(rms, line=inline, xline=xline, 
#                                      time_start=800, time_stop=1500, 
#                                      time_step=2, name='rms')
#            rms_res = rs.resample(trace.time, trace.rms, 
#                                    time_start_re = time_start+c,
#                                    time_stop_re = time_stop-c,
#                                    time_step_re = time_step,
#                                    kind="quadratic", plot=False, 
#                                    logname='rms')        
#            sample = pd.concat([lith_res, 
#                                trace_res.amp,
#                                rms_res.rms], axis=1)
            
            '''generate the dataset for CWT-CNNs'''
            print("==== ==== generating dataset")
            timename = 'time' # name of the time column in the sample  
            yname = 'lith'    # name of the target lith
            xname = ['amp']   # 'rms', 'freq', 'phase', 'IMP'
            
            if i==0:
                x, y, time, seis, freqs = data.cwt2d(
                                       sample, timename, yname, xname, 
                                       num=num, interval=interval,
                                       figsize=figsize, plot=cwtplot)
            else:
                x_w, y_w, time_w, seis_w, freqs_w = data.cwt2d(
                                       sample, timename, yname, xname, 
                                       num=num, interval=interval,
                                       figsize=figsize, plot=cwtplot)
                x = np.append(x, x_w, axis=0)
                y = np.append(y, y_w)
                time = np.append(time, time_w)
                seis = np.append(seis, seis_w)
                #freqs = np.append(freqs, freqs_w)
        '''
        x: cwt map data set
        y: the target parameter
        time: the time of the samples
        seis: the original sesimic values
        freqs: the resampled frequencies
        sample: the original integrated logs and seismic data file
        '''
        return x, y, time, seis, freqs, sample
    
    
    def raw1d(wells, folderpath='data/zj_las', num=10, interval=2):
        '''
        Generating 1d vector data set for DNNs, CNNs
        '''
        '''
        wells: a list of well names for sampling
        folderpath: the folder path for log files
        
        num: number of data beside the targeted data, 
                    total num will be 2*sample+1
        interval: time interval between samples
        '''
        dataset = pd.DataFrame()
        for i in np.arange(0,len(wells)): 
            ## read the las log file
            print("== Generating the 1d dataset of well %s"%wells[i])
            print("==== ==== resampling")
            logs = rlas.las_to_df(file='%s/%s' %(folderpath, wells[i]), 
                                  interval=False,
                                  depth_start=0, depth_stop=0, nan=-999.25)
            
            ## read the time gate and trace location from inf file
            inline = inf.loc[inf.name==wells[i], ['inline']].iloc[0,0]
            xline = inf.loc[inf.name==wells[i], ['xline']].iloc[0,0]
            time_start = inf.loc[inf.name==wells[i], ['start']].iloc[0,0]
            time_stop = inf.loc[inf.name==wells[i], ['stop']].iloc[0,0]
            time_step = 1  #resampling time steps
            
            ## read and resample the log data
            lith_res = rs.resample(logs.TIME, logs.LITH, 
                                   time_start_re = time_start+c,
                                   time_stop_re = time_stop-c,
                                   time_step_re = time_step,
                                   kind="nearest", plot=False, logname='lith')
            
            ## read and resample the seismic data
            trace = rsegy.extract_seismic(segy, line=inline, xline=xline, 
                                      time_start=800, time_stop=1500, 
                                      time_step=2, name='amp')
            trace_res = rs.resample(trace.time, trace.amp, 
                                    time_start_re = time_start+c,
                                    time_stop_re = time_stop-c,
                                    time_step_re = time_step,
                                    kind="quadratic", plot=False, logname='amp')
        
            sample = pd.concat([lith_res, 
                                trace_res.amp], axis=1)
            
#            '''read and resample the rms data'''
#            trace = rsegy.extract_seismic(rms, line=inline, xline=xline, 
#                                      time_start=800, time_stop=1500, 
#                                      time_step=2, name='rms')
#            rms_res = rs.resample(trace.time, trace.rms, 
#                                    time_start_re = time_start+c,
#                                    time_stop_re = time_stop-c,
#                                    time_step_re = time_step,
#                                    kind="quadratic", plot=False, 
#                                    logname='rms')  
#            sample = pd.concat([lith_res, 
#                                trace_res.amp,
#                                rms_res.rms], axis=1)
            
            '''generate the dataset for DNNs, CNNs'''
            print("==== ==== generating dataset")
            time = 'time'    # name of the time column in the sample  
            yname = 'lith'   # name of the target lith
            xname = ['amp']  # 'rms'
            
            data_well = data.raw1d(sample, time, yname, xname,  
                                   num=num, interval=interval)
            dataset = pd.concat([dataset, data_well])
        x = np.array(dataset.iloc[:, 2:((2*num+1)*len(xname)+2)])
        y = np.array(dataset.loc[:, yname]).astype(int)
        return x, y
 
class oversampling:
    '''
    oversampling the dataset
    '''
    def o2d(X, Y, symbol=1, times=2):
        '''
        oversampling the 2d cwt map samples
        '''
        '''
        X: the input maps samples
        Y: the target
        symbol: 1-sand, 0-shale
        times: times to original quantities for oversampling
        '''
        for i in np.arange(0, Y.shape[0]):
            if Y[i] == symbol:
                for j in np.arange(0, (times-1)):
                    Y = np.append(Y, Y[i])
                    X = np.vstack(
                    (X, X[i,:,:,:].reshape(1,X.shape[1],X.shape[2],X.shape[3]))
                                 )
        return X, Y  
    def o1d(X, Y, symbol=1, times=2):
        '''
        oversampling the 1d data set
        '''
        '''
        X: the input samples
        Y: the target
        symbol: 1-sand, 0-shale
        times: times to original quantities for oversampling
        '''
        for i in np.arange(0, Y.shape[0]):
            if Y[i] == symbol:
                for j in np.arange(0, (times-1)):
                    Y = np.append(Y, Y[i])
                    X = np.vstack(
                        (X, X[i,:].reshape(1,X.shape[1]))
                        )
        return X, Y

class line_predict:
    '''
    the whole line processing using established models
    '''
    def predict(model, cwt=True, cnn=True,
                 inline=3662, xline_min=1620, xline_max=1800, 
                 time_start_re = 900, time_stop_re = 1400, time_step_re = 1,
                 num=10, interval=2, channel=1):
        '''
        model: the trainded model
        cwt: True, if using cwt maps
             or False, if using 1d vectors
        cnn: True, if using convoluton layer
             or False, if only using fully connected
        inline: the line number
        xline_min: the min of the xline
        xline_max: the max of the xline
        time_start_re: the start time of the time gate of the line
        time_stop_re: the stop time of the time gate of the line
        time_step_re: the time step resampled to
        num: number of data beside the targeted data, 
                    total num will be 2*sample+1
        interval: time interval between samples
        channel: the channel of the input data
        '''
        for xx in np.arange(xline_min, (xline_max+1)):  #2080
            if xx%50==0: print("==== predicting xline %d trace" %xx)
            trace = rsegy.extract_seismic(segy, line=inline, xline=xx, 
                                      time_start=time_start_re-10, 
                                      time_stop=time_stop_re+10, 
                                      time_step=2, name='amp')
            trace_res = rs.resample(trace.time, trace.amp, 
                                    time_start_re = time_start_re,
                                    time_stop_re = time_stop_re,
                                    time_step_re = time_step_re,
                                    kind="quadratic", plot=False, 
                                    logname='amp')    
            if cwt==True:
                X_seis = data.cwt2d_seis(trace_res, timename='time', 
                            xname=['amp'], num=num, interval=interval, 
                            figsize=(5,5), plot=False)
                X_seis = X_seis.reshape(-1, channel, widths, length)/10000
            else: 
                X_seis = data.raw1d_seis(trace_res, timename='time', 
                                         xname=['amp'],
                                         num=num, interval=interval)
                X_seis = np.array(X_seis)
                if cnn==True: 
                    X_seis = X_seis.reshape(-1, channel, (2*num+1))/10000
                else: X_seis /= 10000
            
            Y_seis = model.predict(X_seis)
            if (xx-xline_min) == 0:
                yline = Y_seis[:,1]
                yline = yline.reshape(-1,1)
            else: 
                yt = Y_seis[:,1].reshape(-1,1)
                yline = np.hstack((yline, yt))
        return yline


    def plot(yline, 
              xline_min=1620, xline_max=1800, 
              time_start_re = 900, time_stop_re = 1400,
              num=10, interval=2, 
              vmax=1, vmin=0, colorbar=True,
              interpolation='nearest', cmap='viridis'):
        '''
        yline:  the predicted result of a whole line
        xline_min: the min of the xline
        xline_max: the max of the xline
        time_start_re: the start time of the time gate of the line
        time_stop_re: the stop time of the time gate of the line
        num: number of data beside the targeted data, 
                    total num will be 2*sample+1
        interval: time interval between samples
        
        vmax: 1 for default
        vmin: 0 for default
        colorbar: True or False
        interpolation: 'quadric', 'nearest', 'bilinear'
        cmap: lithcmap, lithsandshale, 'viridis'
        '''
        plt.figure(figsize=(16,8), facecolor="white")        
        plt.imshow(yline, interpolation=interpolation, aspect='auto', 
                   extent=(xline_min, xline_max, 
                           (time_stop_re-num*interval), 
                           (time_start_re+num*interval)),
                   cmap=cmap, vmax=vmax, vmin=vmin
                   )
        if colorbar==True: plt.colorbar()    


def roc_plot(y_true, y_score, pos_label=None):
    # compute fpr, tpr, thresholds and roc_auc
    fpr, tpr, thresholds = roc_curve(y_true, y_score, pos_label=pos_label)
    roc_auc = auc(fpr, tpr) # compute area under the curve
    plt.figure(figsize=(5,4))
    plt.plot(fpr, tpr, label='ROC curve (AUC area = %0.2f)' % (roc_auc))
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate (FPR)', fontsize=13)
    plt.ylabel('True Positive Rate (TPR)', fontsize=13)
    plt.grid()
    #plt.title('Threshold tuning using ROC')
    plt.legend(loc="lower right")     
    # create the axis of thresholds (scores)
    ax2 = plt.gca().twinx()
    ax2.plot(fpr, thresholds, markeredgecolor='r',linestyle='dashed',color='r')
    ax2.set_ylabel('Threshold',color='r')
    ax2.set_ylim([0,1]) #[thresholds[-1],thresholds[0]]
    ax2.set_xlim([fpr[0],fpr[-1]])
    plt.show()


def to_classes(Y_predict, threshold=0.5):
    '''
    Convert the scores to classes according to the threshold
    '''
    '''
    Y_predicted: the column to convert
    threshold: the threshold to divide the classes
    '''
    Y_predict_type = np.array(Y_predict)
    for i in np.arange(0, Y_predict.shape[0]):
        if Y_predict[i] >= threshold:
            Y_predict_type[i] = 1
        else: Y_predict_type[i] = 0
    return Y_predict_type


class classifier:
    def fnn(X_train, X_test, Y_train, Y_test, activation, neurons, drop, 
            optimizer, batchsize, epoches, classes=2):
        '''
        the deep feedforward neural network model dealing with 1d vector input
        with only fully connected layers
        '''
        model = Sequential()
        model.add(Dense(units=neurons, input_dim=X_train.shape[1]))
        model.add(Activation(activation)) 
        model.add(Dropout(drop))
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))
        model.add(Dense(classes))
        model.add(Activation("softmax"))
        
        #binary_accuracy   categorical_crossentropy
        model.compile(loss = 'binary_crossentropy', 
                      optimizer = optimizer, 
                      metrics= ['accuracy'])
        reduce = callbacks.ReduceLROnPlateau(
                 monitor='val_loss',factor=0.25, patience=10, mode='min', 
                 epsilon=0.005,cooldown=0, min_lr=0.000001)
        earlystop = callbacks.EarlyStopping(monitor='val_loss', 
                                            min_delta=0.001, 
                                            patience=20, 
                                            verbose=0, mode='min')
        print('Training -----------')
        hist = model.fit(X_train, Y_train, 
                         nb_epoch = epoches, 
                         batch_size = batchsize, 
                         validation_data = (X_test, Y_test),
                         callbacks = [reduce, earlystop])
        return model, hist


    def cnn1d(X_train, X_test, Y_train, Y_test, activation, neurons, drop, 
               optimizer, batchsize, epoches, channel=1, classes=2,
               sample_num=10, filters=50, kernel_size=10, strides=2):
        '''
        the convolutional neural networks dealing with 1d vector input
        '''
        model = Sequential()
        model.add(Conv1D(input_shape=(channel, (sample_num*2+1)),
                         filters=filters, kernel_size=kernel_size, 
                         strides=strides, padding='same'))
#        model.add(MaxPooling1D(pool_size=2, strides=2, padding='same'))
#        model.add(Dropout(drop))
        
        model.add(Flatten())
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))
#        model.add(Dense(units=neurons))
#        model.add(Activation(activation))
#        model.add(Dropout(drop))
        
        model.add(Dense(classes))
        model.add(Activation("softmax"))
        
        model.compile(loss = 'binary_crossentropy', 
                      optimizer = optimizer, 
                      metrics= ['accuracy'])
        reduce = callbacks.ReduceLROnPlateau(
                 monitor='val_loss',factor=0.25, patience=10, mode='min', 
                 epsilon=0.005,cooldown=0, min_lr=0.000001)
        earlystop = callbacks.EarlyStopping(monitor='val_loss', 
                                            min_delta=0.001, 
                                            patience=30, 
                                            verbose=0, mode='min')
        print('Training -----------')
        hist = model.fit(X_train, Y_train, 
                         nb_epoch = epoches, 
                         batch_size = batchsize, 
                         validation_data = (X_test, Y_test),
                         callbacks = [reduce, earlystop])
        return model, hist
    
    
    def cnn2d(X_train, X_test, Y_train, Y_test, 
              activation, neurons, drop, optimizer, batchsize, epoches,  
              reduce_threshold=0.01, lr_patience=10, lr_reduce=0.25, 
              stop_threshold=0.01, stop_patience=25,
              filters=25, kernel_size=3, strides=2):
        '''
        the convolutional neural networks dealing with 2d cwt map input
        '''
        model = Sequential()
        model.add(Conv2D(input_shape=X_train[0].shape,
                         filters=filters, kernel_size=kernel_size, 
                         strides=strides, padding='same', 
                         data_format='channels_first'))
        model.add(Flatten())
        
        #model.add(Activation(activation))
#        model.add(MaxPooling2D(pool_size=2, strides=2, padding='same'))
#        model.add(Dropout(drop))
        
#        model.add(Conv2D(filters=filters, kernel_size=kernel_size, 
#                         strides=filters, padding='same',
#                         data_format='channels_first'))
        #model.add(Activation(activation))
        #model.add(MaxPooling2D(pool_size=2, strides=2, padding='same'))
        #model.add(Dropout(drop))      

#        model.add(Conv2D(filters=filters, kernel_size=kernel_size, 
#                         strides=strides, padding='same',
#                         data_format='channels_first'))

        '''used for CWT-DNNs'''
#        model.add(Flatten(input_shape=X_train[0].shape))        
#        model.add(Dense(units=neurons))
#        model.add(Activation(activation))
#        model.add(Dropout(drop))
        
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))       
        
        model.add(Dense(classes))
        model.add(Activation("softmax"))
        
        #categorical_crossentropy ,  binary_crossentropy
        model.compile(loss = 'binary_crossentropy', 
                      optimizer = optimizer, 
                      metrics= ['accuracy'])
        reduce = callbacks.ReduceLROnPlateau(
                monitor='val_loss', factor=lr_reduce, patience=lr_patience, 
                mode='min', 
                epsilon=reduce_threshold, cooldown=0, min_lr=0.000001)
        earlystop = callbacks.EarlyStopping(monitor='val_loss', 
                                            min_delta=stop_threshold, 
                                            patience=stop_patience, 
                                            verbose=0, mode='min')
        print('Training -----------')
        hist = model.fit(X_train, Y_train, 
                         nb_epoch = epoches, 
                         batch_size = batchsize, 
                         validation_data = (X_test, Y_test),
                         callbacks = [reduce, earlystop])
        return model, hist   
    
    def plot(hist):
        '''
        plot the training curve
        
        hist: the training history file
        '''
        cost = pd.DataFrame(hist.history)
        # plotting the loss
        plt.figure(figsize=(12,9), facecolor="white") #dpi=120,
        plt.subplot(2, 2, 1)
        plt.plot(np.arange(1,len(cost.loss)+1), cost.loss,
                 linewidth=1, c='green', label='Training data')
        plt.plot(np.arange(len(cost.val_loss)), cost.val_loss,
                 linewidth=1, c='orange', label='Validating data')
        plt.xlabel("Epoch",fontsize=14)
        plt.ylabel("Loss",fontsize=14)
        plt.legend(loc='upper right', fontsize=12)
        
        # plotting the accuracy
        plt.subplot(2, 2, 2)
        print(cost)
        plt.plot(np.arange(1,len(cost.accuracy)+1), cost.accuracy,
                 linewidth=1, c='green', label='Training data')
        plt.plot(np.arange(len(cost.val_accuracy)), cost.val_accuracy,
                 linewidth=1, c='orange', label='Validating data')
        plt.xlabel("Epoch",fontsize=14)
        plt.ylabel("Accuracy",fontsize=14)
        plt.legend(loc='lower right', fontsize=12)
        
        # plotting the learning rate
        plt.subplot(2, 2, 3)
        plt.plot(np.arange(1,len(cost.lr)+1), cost.lr, 
                 linewidth=1, c='r', label='Learning rate')
        plt.xlabel("Epoch",fontsize=14)
        plt.ylabel("Learning rate",fontsize=14)
        plt.legend(loc='lower left' , fontsize=12)


#%% #################       Main processes    #################################
if __name__ == "__main__":    
    inf = pd.read_excel('data/zj_las/inf.xlsx', 'Sheet1', index_col=None, 
                        na_values=['NA'])    
    print("==== Reading seismic data")
    #segy_high = _read_segy('data/seis/CX_ZJ.SGY', headonly=True)
    segy = _read_segy('data/seis/CX_ZJ_ori.SGY', headonly=True)
    #rms = _read_segy('data/seis/CX_ZJ_RMS.SGY', headonly=True)
    #freq = _read_segy('data/seis/CX_ZJ_Freq.SGY', headonly=True)
    #phase = _read_segy('data/seis/CX_ZJ_Phase.SGY', headonly=True)
    #IMP = _read_segy('data/seis/CX_ZJ_IMP.SGY', headonly=True)
#    segy.traces[0].header
#    rsegy.segy_inf(segy)
#    line3625 = rsegy.extract_line(segy, 3625)
#    rsegy.plot_line(line3625, 
#                    xleft=1620, xright=2080, yup=800, ydown=1500,
#                     cmap="RdBu", percent=99)
    #%%   parameter define  
    '''data processing parameter'''
    wells_train = ['JS8H','ZJ12','LX1', 'CJ186', 'JS2', 'ZJ15',
                   'ZJ11-1', 'ZJ11-3','ZJ13-1', 'ZJ13-3',
                   'JS33-5-1', 'JS33-5-3',  'JS33-5-5',
                   'JS33-2','ZJ10-1', 'ZJ10-3', 'ZJ10-5',
                   'JS33-17-1','JS33-17-3','JS33-17-5',
                   'JS9', 'JS33-15-1', 'JS33-15-3', 'JS33-15-5']
    wells_blind = ['JS7'] #
                   #'JS33-15-4', 'JS33-15-5']
                   #'JS33-17-2', 'JS33-17-4']
                   #'JS33-5-4', 'JS33-5-5']
                   #'ZJ11-3']
                   #'ZJ13-3']
                   #'ZJ10-4', 'ZJ10-5']
    num = 10      # 取样个数  2*num+1
    interval = 2  # interval for sampling 
    widths = 15   # for cwt 10-80hz得到的分频个数
    length = 2*num + 1
    c = 0         # (10-num)*2 ; the time correction for sampling 
    # freqmax=100 # maximum number of frequency for cwt
    
    '''model training parameter'''
    seed = 101        # 101
    testsize = 0.4    # the ratio to split the dataset
    channel=1         # the number of channels 
    classes=2         # the number of classes of the target parameter
    epoches=1000      # the maximum epoches for the training
    
    '''plot parameter'''
    Threshold = 0.4   # the threshold the class the result 
    
    #self-defined color maps.  lithcmap = 'viridis'  'summer' 'Oranges'
    lithcmap = LinearSegmentedColormap.from_list("", 
               ["darkgray",'darkgray','green','greenyellow','yellow','yellow'])
    lithcmap_imp = LinearSegmentedColormap.from_list("", 
                                                     ["yellow","g","darkgray"])
    lithsandshale = ListedColormap(["darkgray","yellow"])   

    #%% 2d cnn -- data  #################################
    '''generate the training and validating data'''
    X, Y, T, seis, freqs, sample = dataset.cwt2d(wells_train, 
                                  folderpath='data/zj_las',
                                  num=num, interval=interval,
                                  figsize=(8,5), cwtplot=False)   
    
    ## split the training and the validation data set
    X_train, X_test, Y_train, Y_test = \
            train_test_split(X, Y, test_size=testsize, random_state=seed)
    
    ## oversample the sands
    X_train, Y_train = oversampling.o2d(X_train, Y_train, symbol=1, times=2)
    X_test, Y_test = oversampling.o2d(X_test, Y_test, symbol=1, times=2)    
    
    ## data preprocessing    
    Y_train = np_utils.to_categorical(Y_train, num_classes=classes)
    Y_test = np_utils.to_categorical(Y_test, num_classes=classes)
    X_train = X_train.reshape(-1, channel, widths, length)/10000
    X_test = X_test.reshape(-1, channel, widths, length)/10000

    '''generate the blind test data'''
    X_blind, Y_blind, T_blind, seis_blind, freqs_b, sample_b = dataset.cwt2d(
                                            wells_blind,
                                            folderpath='data/zj_las',
                                            num=num, interval=interval,
                                            figsize=(8,5), cwtplot=False)
    Y_blind = np_utils.to_categorical(Y_blind, num_classes=classes)
    X_blind = X_blind.reshape(-1, channel, widths, length) /10000   
    #%% 2d cnn -- train 
    model, hist = classifier.cnn2d(X_train, X_test, Y_train, Y_test,  
                       activation='relu', neurons=80, drop=0.6, 
                       optimizer=Adam(lr=0.0001),  #'SGD' Adam(lr=0.0001)
                       batchsize=40, epoches=epoches,
                       reduce_threshold=0.005, lr_patience=10, lr_reduce=0.25,
                       stop_threshold=0.005, stop_patience=25,
                       filters=40, kernel_size=7, strides=2)
    classifier.plot(hist)

   #save the trained model
#    model.save('model/model_fnn-cwt_20200131.h5')
    model.save('model/model_cnn2d_20200131_cwt-cnn.h5')
    cost = pd.DataFrame(hist.history)
    cost.to_csv('model/history_cnn2d_20200131_cwt-cnn.csv', 
               index_label='epoch')

#    #load the saved model
#    model = load_model('model/model_fnn-cwt_20200131.h5')
#    from keras.utils import plot_model
#    plot_model(model, to_file='model/model_cnn2d_20200131-good.png')

    #%% 2d cnn -- analysis    
    layer_list = list([(layer.name, layer) for layer in model.layers])    
    
    #model.summary()

    # the interlayer output of CWT-CNNs
    layer_name = layer_list[2][0]  #'block5_conv3' 'conv2d_1'
    inter_layer_model = Model(inputs=model.input,
                              outputs=model.get_layer(layer_name).output)
    inter_output = inter_layer_model.predict(X_blind)        
    #%%  #################################
    '''plot the input maps '''
    freq = [81.4,74.0,67.8,62.6,58.1,54.2,50.8,45.2,40.7,35.4,
            30.1,25.4,20.3,15.0,10.0]
    
    start=960    # the start time of the sampled well dataset
    target=1170  # the target time of the sampled well dataset to plot
    
    cwtmatr = X_blind[target-start,0,:,:]  #start=960, target=1074
    cwtmatrT = cwtmatr.T
    maximum = 1
    cwtmatrT = np.clip(cwtmatrT, -maximum, maximum, cwtmatrT)
    levels = np.linspace(-maximum, maximum, 41)
    plt.figure(figsize=(5,4))
    plt.contourf(freq, np.arange(target-20, target+21, 2),  cwtmatrT, 
                 levels=levels, cmap='RdBu')
    #plt.ylim(10, 80)
    ax = plt.gca()
    ax.invert_yaxis()
    plt.xlabel('Time(ms)')
    plt.ylabel('Frequency(Hz)')
    plt.colorbar()
    
    #%%  #################################
    '''plot interlayer feature maps'''
    q = 20 #5 10 15 20
    plt.figure(figsize=(10,2*q/5))
    """
    for i in np.arange(0,q,1):
        plt.subplot(q/5,5,i+1)
        cwtmatr2 = inter_output[target-start,i,:,:] 
        maximum = 1
        np.clip(cwtmatr2, -maximum, maximum, cwtmatr2)
        levels = np.linspace(-maximum, maximum, 41)        
        plt.contourf(np.arange(0, inter_output.shape[2]), 
                     np.arange(0, inter_output.shape[3]), 
                     cwtmatr2.T, 
                     levels=levels, cmap='RdBu')
        #plt.ylim(10, 80)
        ax = plt.gca()
        ax.invert_yaxis()
        ax.invert_xaxis()
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_xticks([])
        plt.tight_layout()
        #plt.colorbar()
    """
    #%% 2d cnn -- threshold
    '''testing on blind data''' 
    loss, accuracy = model.evaluate(X_blind, Y_blind)
    print('\ntest loss: ', loss)
    print('test accuracy: ', accuracy)
    Y_predict = model.predict(X_blind)
    Y_train_predict = model.predict(X_train)
    Y_test_predict = model.predict(X_test)
    
    '''Threshold tuning using ROC'''
    roc_plot(Y_train[:,1], Y_train_predict[:,1], pos_label=1)
    roc_plot(Y_test[:,1], Y_test_predict[:,1], pos_label=1)
    roc_plot(Y_blind[:,1], Y_predict[:,1], pos_label=1)
    
    '''Threshold tuning using F1 curve''' 
    F1_blind_c = []
    F1_train_c = []
    F1_test_c = []
    for thr in np.arange(0,1,0.01):
        Y_pred_type = to_classes(Y_predict[:,1], threshold=thr)
        Y_train_pred_type = to_classes(Y_train_predict[:,1], threshold=thr)
        Y_test_pred_type = to_classes(Y_test_predict[:,1], threshold=thr)
        F1_blind = f1_score(Y_blind[:,1], Y_pred_type, average='weighted') 
        F1_train = f1_score(Y_train[:,1], Y_train_pred_type,average='weighted')
        F1_test = f1_score(Y_test[:,1], Y_test_pred_type, average='weighted')
        F1_blind_c.append(F1_blind)
        F1_train_c.append(F1_train)
        F1_test_c.append(F1_test)
    
    plt.figure(figsize=(7,5))
    plt.plot(np.arange(0,1,0.01), F1_train_c, 'b', label='Training data')
    plt.plot(np.arange(0,1,0.01), F1_test_c, 'g', label='Validating data')
    plt.plot(np.arange(0,1,0.01), F1_blind_c, 'r', label='Testing data')
    plt.xlim(0,1)
    plt.xlabel('Threshold', fontsize=13)
    plt.ylabel('F1 score', fontsize=13)
    plt.legend(loc='lower right')
    plt.grid()
    
#    np.argmax(F1_blind_c)
#    np.argmax(F1_test_c)
#    np.argmax(F1_train_c)
    #%% 2d cnn -- plot on blind data
    threshold=0.4
    Y_pred_type = to_classes(Y_predict[:,1], threshold=threshold)
    Y_train_pred_type = to_classes(Y_train_predict[:,1], threshold=threshold)
    Y_test_pred_type = to_classes(Y_test_predict[:,1], threshold=threshold)
    print("Classification report for the blind data:\n%s\n"
          % metrics.classification_report(Y_blind[:,1], Y_pred_type))
    print("Confusion matrix for the blind data:\n%s" 
          % metrics.confusion_matrix(Y_blind[:,1], Y_pred_type))   
    print("Classification report for the training data:\n%s\n"
          % metrics.classification_report(Y_train[:,1], Y_train_pred_type))
    print("Confusion matrix for the training data:\n%s" 
          % metrics.confusion_matrix(Y_train[:,1], Y_train_pred_type))        
    print("Classification report for the validating data:\n%s\n"
          % metrics.classification_report(Y_test[:,1], Y_test_pred_type))
    print("Confusion matrix for the validating data:\n%s" 
          % metrics.confusion_matrix(Y_test[:,1], Y_test_pred_type))        
    print("\nF1 for the blind, validation, training data: %.3f %.3f %.3f"
          % (f1_score(Y_blind[:,1], Y_pred_type, average='weighted'), 
            f1_score(Y_test[:,1], Y_test_pred_type, average='weighted'),
            f1_score(Y_train[:,1], Y_train_pred_type, average='weighted'))
          ) 
          
    Y_target = np.stack((Y_blind[:,1], Y_blind[:,1]), axis=-1)
    Y_comp = np.stack((Y_pred_type, Y_predict[:,1]), axis=-1)    
    
    tmin=int(T_blind.min(axis=0))
    tmax=int(T_blind.max(axis=0))
    
    plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(1, 5, width_ratios=[1.5, 1, 1.6, 1.1, 2])
    plt.subplot(gs[0])
    plt.plot(seis_blind/10000, T_blind, 'k')

    """
    plt.fill_betweenx(T_blind.iloc[:,0], seis_blind.iloc[:,0]/10000, x2=0, 
                      where=seis_blind.iloc[:,0]/10000 >= 0, 
                      facecolor='g', interpolate=True)
    """               
                      
#    '''for JS33-7'''
#    plt.fill_betweenx(T_blind, seis_blind/10000, x2=0, 
#                      where=seis_blind/10000 >= 0, 
#                      facecolor='k', interpolate=True)

#    T_b_train, T_b_test, Y_b_train, Y_b_test = \
#      train_test_split(T_blind, Y_blind, test_size=testsize, random_state=seed)
#    plt.scatter(Y_b_train[:,1], T_b_train)
    
#    seis_blind2 = seis_blind+20000
#    seis_blind3 = seis_blind-20000    
#    plt.plot(seis_blind2, T_blind, 'k')
#    plt.fill_betweenx(T_blind.iloc[:,0], seis_blind2.iloc[:,0], x2=20000, 
#                      where=seis_blind2.iloc[:,0] >= 20000, facecolor='k', 
#                      interpolate=True)
#    plt.plot(seis_blind3, T_blind, 'k')
#    plt.fill_betweenx(T_blind.iloc[:,0], seis_blind3.iloc[:,0], x2=-20000, 
#                      where=seis_blind3.iloc[:,0] >= -20000, facecolor='k', 
#                      interpolate=True)
    plt.xlim(xmin=-5,xmax=5) 
    plt.ylim(ymin=tmin,ymax=tmax)
    plt.ylabel('Time (ms)', fontsize=13)
    plt.xlabel('Seismic', fontsize=13)
    #plt.yticks(np.arange(tmin,tmax,50))
    plt.subplot(gs[0]).invert_yaxis()

    
    plt.subplot(gs[1])
    plt.imshow(Y_target, interpolation='nearest', aspect='auto', 
               extent=(0, 1, tmax, tmin),
               cmap=lithsandshale)
    plt.xticks([0.5], ['Lith'], fontsize=11)
    plt.xlabel('Target lithology', fontsize=13)    
    #plt.colorbar()
    plt.subplot(gs[1]).set_yticklabels([])
    


    plt.subplot(gs[2])
    plt.imshow(Y_comp, interpolation='nearest', aspect='auto', 
               extent=(0, 2, tmax, tmin), vmin=0,vmax=1,
               cmap=lithcmap)
    plt.xticks([0.5,1.5], ['Lith','Scores'], fontsize=11) 
    plt.xlabel('CNN result', fontsize=13)
    #plt.colorbar()
    plt.subplot(gs[2]).set_yticklabels([])
    
    #%%
#    #do fnn first
#    plt.subplot(gs[3])
#    plt.imshow(Y1d_comp, interpolation='nearest', aspect='auto', 
#               extent=(0, 2, tmax, tmin),
#               cmap=lithcmap)
#    plt.xticks([0.5,1.5], ['Lith','Scores'], fontsize=11) 
#    plt.xlabel('FNN result', fontsize=13)
#    plt.colorbar()
#    plt.subplot(gs[3]).set_yticklabels([])
#
#    plt.subplot(gs[4])
#    trace_imp = rsegy.extract_seismic(IMP, line=3662, xline=1938, 
#                              time_start=800, time_stop=1500, 
#                              time_step=2, name='imp')
#    trace_imp_res = rs.resample(trace_imp.time, trace_imp.imp, 
#                            time_start_re = 963,
#                            time_stop_re = 1358,
#                            time_step_re = 1,
#                            kind="quadratic", plot=False, 
#                            logname='imp')
#    #for JS7: line=3662, xline=1938,time_start_re = 963,time_stop_re = 1358
#    #for JS9: line=3715, xline=1840,time_start_re = 1003,time_stop_re = 1263
#    
#    '''#find the threshold for p-impedence
#    trace_imp_res_type == trace_imp_res
#    F1_imp_c = []
#    for thr in np.arange(8000,12000,100):
#        trace_imp_res_type.loc[trace_imp_res.imp <= thr, 'imp'] = 1
#        trace_imp_res_type.loc[trace_imp_res.imp > thr, 'imp'] = 0
#        F1_imp = f1_score(Y_blind[:,1], trace_imp_res_type.imp,average='weighted')
#        F1_imp_c.append(F1_imp)
#    plt.plot(np.arange(8000,12000,100), F1_imp_c, label='P-impedence')
#    plt.xlabel('Threshold', fontsize=13)
#    plt.ylabel('F1 score', fontsize=13)
#    plt.legend(loc='lower left')
#    plt.grid()  
#    np.argmax(F1_imp_c)
#    '''
#    trace_imp_res_type = trace_imp_res
#    trace_imp_res_type.loc[trace_imp_res.imp <= 9750, 'imp'] = 0
#    trace_imp_res_type.loc[trace_imp_res.imp > 9750, 'imp'] = 100000
#    imp_comp = pd.concat([trace_imp_res_type.imp, trace_imp_res.imp], axis=1)   
#    plt.imshow(imp_comp, interpolation='nearest', aspect='auto', 
#               extent=(0, 2, tmax, tmin), vmin=8500, vmax=11000,
#               cmap=lithcmap_imp)
#    plt.xticks([0.5,1.5], ['Lith','P-impedence'], fontsize=11) 
#    plt.xlabel('Seismic inversion', fontsize=13)
#    plt.colorbar()
#    plt.subplot(gs[4]).set_yticklabels([])
#    # for JS9: 10000, 9000, 11500
    
#    plt.subplot(gs[3])
#    JS7 = rlas.las_to_df(file='data/zj_las/JS33-17-1', interval=False,
#                      depth_start=0, depth_stop=0, nan=-999.25)
#    plt.plot(JS7.GR, JS7.TIME, 'r', linewidth=0.75)
#    plt.xlim(xmin=40,xmax=130) 
#    plt.ylim(ymin=tmin,ymax=tmax)
#    #plt.xlabel('GR', fontsize=13)
#    plt.subplot(gs[3]).invert_yaxis()
#    plt.subplot(gs[3]).set_yticklabels([])
#    
#    plt.tight_layout() 
#    plt.show()
    
    #%% 2d cnn -- hyper parameters
    Y_predict_all = []
    Y_train_predict_all = []
    Y_test_predict_all = []
    F1_blind_all = []
    F1_train_all = []
    F1_test_all = []
    thr_range = np.arange(0,1,0.01)
    for t in np.arange(0,5):
        model, hist = classifier.cnn2d(X_train, X_test, Y_train, Y_test,  
                       activation='relu', neurons=80, drop=0.5, 
                       optimizer=Adam(lr=0.0001),  #'SGD' Adam(lr=0.0001)
                       batchsize=40, epoches=epoches,
                       filters=60, kernel_size=5, strides=2)
        classifier.plot(hist)
        
        Y_predict = model.predict(X_blind)
        Y_train_predict = model.predict(X_train)
        Y_test_predict = model.predict(X_test)
        
        F1_blind_c = []
        F1_train_c = []
        F1_test_c = []
        for thr in thr_range:
            Y_pred_type = to_classes(Y_predict[:,1], threshold=thr)
            Y_train_pred_type = to_classes(Y_train_predict[:,1], threshold=thr)
            Y_test_pred_type = to_classes(Y_test_predict[:,1], threshold=thr)
            F1_blind = f1_score(Y_blind[:,1], Y_pred_type, 
                                average='weighted') #
            F1_train = f1_score(Y_train[:,1], Y_train_pred_type, 
                                average='weighted')
            F1_test = f1_score(Y_test[:,1], Y_test_pred_type, 
                               average='weighted')
            F1_blind_c.append(F1_blind)
            F1_train_c.append(F1_train)
            F1_test_c.append(F1_test)
        if t==0:
            Y_predict_all = Y_predict[:,1]
            Y_train_predict_all = Y_train_predict[:,1]
            Y_test_predict_all = Y_test_predict[:,1]
            F1_blind_all = F1_blind_c
            F1_train_all = F1_train_c
            F1_test_all = F1_test_c
        else:
            Y_predict_all = np.column_stack((Y_predict_all,Y_predict[:,1]))
            Y_train_predict_all = np.column_stack((Y_train_predict_all,
                                                   Y_train_predict[:,1]))
            Y_test_predict_all = np.column_stack((Y_test_predict_all,
                                                  Y_test_predict[:,1]))
            F1_blind_all = np.column_stack((F1_blind_all,F1_blind_c))
            F1_train_all = np.column_stack((F1_train_all,F1_train_c))
            F1_test_all = np.column_stack((F1_test_all,F1_test_c))
    
    #%%
    '''plot the f1 scores'''
    Y_predict_all_mean = np.mean(Y_predict_all, axis=1)
    Y_train_predict_all_mean = np.mean(Y_train_predict_all, axis=1)
    Y_test_predict_all_mean = np.mean(Y_test_predict_all, axis=1)
    Y_predict_all_std = np.std(Y_predict_all, axis=1)
    Y_train_predict_all_std = np.std(Y_train_predict_all, axis=1)
    Y_test_predict_all_std = np.std(Y_test_predict_all, axis=1)
    
    F1_blind_all_mean = np.mean(F1_blind_all, axis=1)
    F1_blind_all_std = np.std(F1_blind_all, axis=1)
    F1_train_all_mean = np.mean(F1_train_all, axis=1)
    F1_train_all_std = np.std(F1_train_all, axis=1)    
    F1_test_all_mean = np.mean(F1_test_all, axis=1)
    F1_test_all_std = np.std(F1_test_all, axis=1)
    
    plt.figure(figsize=(8,5))
    plt.plot(thr_range, F1_train_all_mean, 'b', label='Training data')
    plt.fill_between(thr_range, F1_train_all_mean - F1_train_all_std,
                 F1_train_all_mean + F1_train_all_std, alpha=0.2, color="b")
    plt.plot(thr_range, F1_test_all_mean, 'g', label='Validating data')
    plt.fill_between(thr_range, F1_test_all_mean - F1_test_all_std,
                 F1_test_all_mean + F1_test_all_std, alpha=0.2, color="g")
#    plt.plot(thr_range, F1_blind_all_mean, 'r', label='Testing data')
#    plt.fill_between(thr_range, F1_blind_all_mean - F1_blind_all_std,
#                 F1_blind_all_mean + F1_blind_all_std, alpha=0.2, color="r")
    #plt.xlim(0.0, 0.8)
    #plt.ylim(0.4, 1.0)
    plt.xlabel('Threshold', fontsize=13)
    plt.ylabel('F1', fontsize=13)
    plt.legend(loc='lower right')
    plt.grid()
    plt.show()
    
    plt.figure(figsize=(15,2))
    plt.plot(T_blind, Y_blind[:,1], 'r', label='Training data')
    plt.plot(T_blind, Y_predict_all_mean, 'b', label='Training data')
    plt.fill_between(T_blind.time, Y_predict_all_mean - Y_predict_all_std,
                 Y_predict_all_mean + Y_predict_all_std, alpha=0.25, color="b")
    plt.show()
    
    #%%
    '''plot the result'''
    threshold=0.4
    Y_pred_type = to_classes(Y_predict_all_mean, threshold=threshold)
    Y_train_pred_type = to_classes(Y_train_predict_all_mean, threshold=threshold)
    Y_test_pred_type = to_classes(Y_test_predict_all_mean, threshold=threshold)
    print("Classification report for the blind data:\n%s\n"
          % metrics.classification_report(Y_blind[:,1], Y_pred_type))
    print("Confusion matrix for the blind data:\n%s" 
          % metrics.confusion_matrix(Y_blind[:,1], Y_pred_type))   
    print("Classification report for the training data:\n%s\n"
          % metrics.classification_report(Y_train[:,1], Y_train_pred_type))
    print("Confusion matrix for the training data:\n%s" 
          % metrics.confusion_matrix(Y_train[:,1], Y_train_pred_type))        
    print("Classification report for the validating data:\n%s\n"
          % metrics.classification_report(Y_test[:,1], Y_test_pred_type))
    print("Confusion matrix for the validating data:\n%s" 
          % metrics.confusion_matrix(Y_test[:,1], Y_test_pred_type))        
    print("\nF1 for the blind, validation, training data: %.3f %.3f %.3f"
          % (f1_score(Y_blind[:,1], Y_pred_type, average='weighted'), 
            f1_score(Y_test[:,1], Y_test_pred_type, average='weighted'),
            f1_score(Y_train[:,1], Y_train_pred_type, average='weighted'))
          ) 
    
    Y_target = np.stack((Y_blind[:,1], Y_blind[:,1]), axis=-1)
    Y_comp = np.stack((Y_pred_type, Y_predict[:,1]), axis=-1)    
    
    tmin=int(T_blind.min(axis=0))
    tmax=int(T_blind.max(axis=0))
    
    plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.5, 1, 2])
    plt.subplot(gs[0])
    plt.plot(seis_blind/10000, T_blind, 'k')
    plt.xlim(xmin=-5,xmax=5) 
    plt.ylim(ymin=tmin,ymax=tmax)
    plt.ylabel('Time (ms)', fontsize=13)
    plt.xlabel('Seismic', fontsize=13)
    plt.subplot(gs[0]).invert_yaxis()
   
    plt.subplot(gs[1])
    plt.imshow(Y_target, interpolation='nearest', aspect='auto', 
               extent=(0, 1, tmax, tmin),
               cmap=lithsandshale)
    plt.xticks([0.5], ['Lith'], fontsize=11)
    plt.xlabel('Target lithology', fontsize=13)    
    #plt.colorbar()
    plt.subplot(gs[1]).set_yticklabels([])

    plt.subplot(gs[2])
    plt.imshow(Y_comp, interpolation='nearest', aspect='auto', 
               extent=(0, 2, tmax, tmin),
               cmap=lithcmap)
    plt.xticks([0.5,1.5], ['Lith','Scores'], fontsize=11) 
    plt.xlabel('CNN result', fontsize=13)
    plt.colorbar()
    plt.subplot(gs[2]).set_yticklabels([])
    plt.tight_layout() 
    plt.show()


    #%%  1d cnn -- data  #################################
    '''generate the training and validating data'''
    X1d, Y1d = dataset.raw1d(wells_train, num=num, interval=interval) 
    X1d_train1, X1d_test1, Y1d_train1, Y1d_test1 = \
        train_test_split(X1d, Y1d, test_size=testsize, random_state = seed)    

    X1d_train, Y1d_train = oversampling.o1d(X1d_train1, Y1d_train1, symbol=1, 
                                            times=2) 
    X1d_test, Y1d_test = oversampling.o1d(X1d_test1, Y1d_test1, symbol=1, 
                                          times=2)
    
    X1d_train = X1d_train.reshape(-1, channel, length)/10000
    X1d_test = X1d_test.reshape(-1, channel, length)/10000    
    Y1d_train = np_utils.to_categorical(Y1d_train, num_classes=classes)
    Y1d_test = np_utils.to_categorical(Y1d_test, num_classes=classes)

    '''generate the blind testing data'''
    X1d_blind, Y1d_blind = dataset.raw1d(wells_blind,num=num,interval=interval)
    X1d_blind = X1d_blind.reshape(-1, channel, length)/10000
    Y1d_blind = np_utils.to_categorical(Y1d_blind, num_classes=classes)
    #%% 1d cnn -- train
    model1d, hist1d = classifier.cnn1d(X1d_train, X1d_test,
                                           Y1d_train, Y1d_test,
                       activation='relu', neurons=80, 
                       drop=0.5, optimizer=Adam(lr=0.0001),  # 'SGD'
                       batchsize=40, epoches=epoches,
                       channel=channel,  #the number of x para
                       classes=classes,   #the number of lithologies
                       sample_num=num,
                       filters=60, kernel_size=5, strides=2)

    classifier.plot(hist1d)
    #model1d.save('model/model_cnn1d_20180209.h5')
    #%%  1d fnn -- data   #################################
    '''generate the training and validating data'''
    X1d, Y1d = dataset.raw1d(wells_train, num=num, interval=interval)
    X1d_train, X1d_test, Y1d_train, Y1d_test = \
        train_test_split(X1d, Y1d, test_size=testsize, random_state = seed)

    X1d_train, Y1d_train = oversampling.o1d(X1d_train1, Y1d_train1, 
                                            symbol=1, times=2) 
    X1d_test, Y1d_test = oversampling.o1d(X1d_test1, Y1d_test1, 
                                            symbol=1, times=2)
             
    X1d_train = X1d_train/10000
    X1d_test = X1d_test/10000   
    Y1d_train = np_utils.to_categorical(Y1d_train, num_classes=classes)
    Y1d_test = np_utils.to_categorical(Y1d_test, num_classes=classes)
    
    '''generate the blind testing data'''
    X1d_blind, Y1d_blind = dataset.raw1d(wells_blind,num=num,interval=interval)
    X1d_blind = X1d_blind/10000
    Y1d_blind = np_utils.to_categorical(Y1d_blind, num_classes=classes)
    #%% 1d fnn -- train
    model1d, hist1d = classifier.fnn(X1d_train,X1d_test,Y1d_train,Y1d_test,
                       activation='relu', neurons=80, 
                       drop=0.4, 
                       optimizer=Adam(lr=0.0001),  #'SGD' Adam(lr=0.0001)
                       batchsize=40, 
                       epoches=epoches,
                       classes=classes    #the number of lithologies
                       )
    classifier.plot(hist1d)
    #model1d.save('model/model_fnn_20180210.h5')
       
    #%%  1d fnn and cnn -- threshold
    # model1d = load_model('model/model_cnn1d_20180209.h5')
    loss, accuracy = model1d.evaluate(X1d_blind, Y1d_blind)
    print('\ntest loss: ', loss)
    print('\ntest accuracy: ', accuracy)
    Y1d_predict = model1d.predict(X1d_blind)
    Y1d_train_predict = model1d.predict(X1d_train)
    Y1d_test_predict = model1d.predict(X1d_test)    
        
    '''Threshold tuning using ROC'''
    roc_plot(Y1d_train[:,1], Y1d_train_predict[:,1], pos_label=1)
    roc_plot(Y1d_test[:,1], Y1d_test_predict[:,1], pos_label=1)
    roc_plot(Y1d_blind[:,1], Y1d_predict[:,1], pos_label=1)
    
    '''ploting the results''' 
    F11d_blind_c = []
    F11d_train_c = []
    F11d_test_c = []
    for thr in np.arange(0,1,0.01):
        Y1d_pred_type = to_classes(Y1d_predict[:,1], threshold=thr)
        Y1d_train_pred_type = to_classes(Y1d_train_predict[:,1], threshold=thr)
        Y1d_test_pred_type = to_classes(Y1d_test_predict[:,1], threshold=thr)
        F11d_blind = f1_score(Y1d_blind[:,1], Y1d_pred_type, 
                              average='weighted') #
        F11d_train = f1_score(Y1d_train[:,1], Y1d_train_pred_type, 
                              average='weighted')
        F11d_test = f1_score(Y1d_test[:,1], Y1d_test_pred_type, 
                             average='weighted')
        F11d_blind_c.append(F11d_blind)
        F11d_train_c.append(F11d_train)
        F11d_test_c.append(F11d_test)

    plt.plot(np.arange(0,1,0.01), F11d_blind_c, label='Testing data')
    plt.plot(np.arange(0,1,0.01), F11d_test_c, label='Validating data')
    plt.plot(np.arange(0,1,0.01), F11d_train_c, label='Training data')
    plt.xlabel('Threshold', fontsize=13)
    plt.ylabel('F1 score', fontsize=13)
    plt.legend(loc='lower right')
    plt.grid()
    
#    np.argmax(F1_blind_c)
#    np.argmax(F1_test_c)
#    np.argmax(F1_train_c)
    #%%  1d fnn and cnn -- plot on blind data
    threshold1d=0.4
    Y1d_pred_type = to_classes(Y1d_predict[:,1], threshold=threshold1d)
    Y1d_train_pred_type = to_classes(Y1d_train_predict[:,1], 
                                     threshold=threshold1d)
    Y1d_test_pred_type = to_classes(Y1d_test_predict[:,1], 
                                    threshold=threshold1d)
    print("Classification report for the 1d blind data:\n%s\n"
          % metrics.classification_report(Y1d_blind[:,1], Y1d_pred_type))
    print("Confusion matrix for the 1d blind data:\n%s" 
          % metrics.confusion_matrix(Y1d_blind[:,1], Y1d_pred_type))   
    print("Classification report for the 1d training data:\n%s\n"
          % metrics.classification_report(Y1d_train[:,1], Y1d_train_pred_type))
    print("Confusion matrix for the 1d training data:\n%s" 
          % metrics.confusion_matrix(Y1d_train[:,1], Y1d_train_pred_type))        
    print("Classification report for the 1d validating data:\n%s\n"
          % metrics.classification_report(Y1d_test[:,1], Y1d_test_pred_type))
    print("Confusion matrix for the 1d validating data:\n%s" 
          % metrics.confusion_matrix(Y1d_test[:,1], Y1d_test_pred_type))    
    print("\nF1 for the blind, validation, training data: %.3f %.3f %.3f"
          % (f1_score(Y1d_blind[:,1], Y1d_pred_type, average='weighted'), 
            f1_score(Y1d_test[:,1], Y1d_test_pred_type, average='weighted'),
            f1_score(Y1d_train[:,1], Y1d_train_pred_type, average='weighted'))
          ) 

    Y_target = np.stack((Y_blind[:,1], Y_blind[:,1]), axis=-1)
    Y1d_comp = np.stack((Y1d_pred_type, Y1d_predict[:,1]), axis=-1)    
    '''ploting the results'''
    tmin=int(T_blind.min(axis=0))
    tmax=int(T_blind.max(axis=0))
    
    plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(1, 5, width_ratios=[1.5, 1, 1.6, 2, 2])
    plt.subplot(gs[0])
    plt.plot(seis_blind/10000, T_blind, 'k')
#    plt.fill_betweenx(T_blind.iloc[:,0], seis_blind.iloc[:,0], x2=0, 
#                      where=seis_blind.iloc[:,0] >= 0, facecolor='k', 
#                      interpolate=True)

#    seis_blind2 = seis_blind+20000
#    seis_blind3 = seis_blind-20000    
#    plt.plot(seis_blind2, T_blind, 'k')
#    plt.fill_betweenx(T_blind.iloc[:,0], seis_blind2.iloc[:,0], x2=20000, 
#                      where=seis_blind2.iloc[:,0] >= 20000, facecolor='k', 
#                      interpolate=True)
#    plt.plot(seis_blind3, T_blind, 'k')
#    plt.fill_betweenx(T_blind.iloc[:,0], seis_blind3.iloc[:,0], x2=-20000, 
#                      where=seis_blind3.iloc[:,0] >= -20000, facecolor='k', 
#                      interpolate=True)
    
    plt.xlim(xmin=-5,xmax=5)
    plt.ylim(ymin=tmin,ymax=tmax)
    plt.ylabel('Time (ms)', fontsize=13)
    plt.xlabel('Seismic', fontsize=13)
    plt.subplot(gs[0]).invert_yaxis()

    plt.subplot(gs[1])
    plt.imshow(Y_target, interpolation='nearest', aspect='auto', 
               extent=(0, 1, tmax, tmin),
               cmap=lithsandshale)
    plt.xticks([0.5], ['Lith'], fontsize=11)
    plt.xlabel('Target lithology', fontsize=13)    
    #plt.colorbar()
    plt.subplot(gs[1]).set_yticklabels([])

    plt.subplot(gs[2])
    plt.imshow(Y1d_comp, interpolation='nearest', aspect='auto', 
               cmap=lithcmap, vmin=0,vmax=1,
               extent=(0, 2, tmax, tmin)) #'summer'
    plt.xticks([0.5,1.5], ['Lith','Scores'], fontsize=11) 
    plt.xlabel('FNN result', fontsize=13)
    #plt.colorbar()
    plt.subplot(gs[2]).set_yticklabels([])
     
    #%% apply model to a whole line   #######################
    model = load_model('model/model_cnn2d_20200131_cwt-cnn.h5')
    Yline = line_predict.predict(model, cwt=True, cnn=True,
                 inline=3662, xline_min=1620, xline_max=2080,  #2080
                 time_start_re = 900, time_stop_re = 1400, time_step_re = 1,
                 num=num, interval=interval, channel=channel)
    pd.DataFrame(Yline).to_csv("output/Yline_fnn-cwt_20200131.csv")

    #%% plot the imported data
    Ylinei = pd.read_csv('output/Yline_fnn-cwt_20200131.csv', index_col=None, 
                        na_values=['NA']) 
    Ylinei = np.array(Ylinei)
    Ylinei = np.delete(Ylinei, 0, 1)
    Ylinei[Ylinei >= 0.42] = 1
    Ylinei[Ylinei < 0.42] = 0
    
    line_predict.plot(Ylinei, 
                  xline_min=1620, xline_max=2080, 
                  time_start_re = 900, time_stop_re = 1400,
                  num=num, interval=interval,
                  interpolation='bilinear', 
                  cmap=lithsandshale,  #lithcmap   lithsandshale
                  colorbar=False)

    #%% plot the seismic profile
    line2 = rsegy.extract_line_v2(segy, 
                                inline=3662, xline_min=1620, xline_max=2080, 
                                time_start=900, time_stop=1400, time_step=2)

    line_predict.plot(line2/10000, 
                  xline_min=1620, xline_max=2080, 
                  time_start_re = 900, time_stop_re = 1400,
                  num=num, interval=interval, vmax=2, vmin=-2,
                  interpolation='bilinear', cmap='RdBu', colorbar=True)
    #'viridis'
    #%% apply model to a 3D cube      #######################

    model = load_model('model/model_cnn2d_20200131_cwt-cnn.h5')
    
    time_start = 900  #900
    time_stop = 1400  #1400
    xline_min=1620 
    xline_max=2080
    inline_min = 3120  #3120  JS7-3662
    inline_max = 3800  #3800
    Ycube = np.zeros((inline_max-inline_min,
                      time_stop-time_start-num*interval*2,
                      xline_max-xline_min+1))

    for i in np.arange(inline_min, inline_max):
        if i%2==0: print("==== predicting inline %d traces ====" %i)
        Yline = line_predict.predict(model, cwt=True, cnn=True,
                 inline=i, xline_min=xline_min, xline_max=xline_max,  #2080
                 time_start_re = time_start, time_stop_re = time_stop, 
                 time_step_re = 1,num=num, interval=interval, channel=channel)
        Ycube[i-inline_min,:,:] = Yline
    #export
    Ycube1 = Ycube.reshape(Ycube.shape[0],Ycube.shape[1]*Ycube.shape[2])
    np.savetxt("output/Ycube_cnn-cwt_20200131.csv", Ycube1, delimiter=',')
    #import
    Ycube2 = pd.read_csv("output/Ycube_cnn-cwt_20200131.csv", index_col=None, 
                         header=None, na_values=['NA'])
    Ycube2 = np.array(Ycube2).reshape(inline_max-inline_min,
                                      time_stop-time_start-num*interval*2,
                                      xline_max-xline_min+1)
    
    
    '''plot the verticle slice'''
    line_predict.plot(Ycube2[0,:,:], 
                  xline_min=xline_min, xline_max=xline_max, 
                  time_start_re = time_start, time_stop_re = time_stop,
                  num=num, interval=interval,
                  interpolation='bilinear', 
                  cmap=lithcmap,  #lithcmap   lithsandshale
                  colorbar=False)

    '''plot the horizontal slice'''
    plt.figure(figsize=(8,8), facecolor="white")        
    plt.imshow(Ycube2[:,5,:], interpolation='bilinear', aspect='auto', 
               extent=(xline_min, xline_max, inline_min, inline_max),
               cmap=lithcmap, vmax=1, vmin=0
               )

    












    
  