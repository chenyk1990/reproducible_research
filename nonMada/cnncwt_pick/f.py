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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import struct
 
from sklearn.metrics import accuracy_score, f1_score
from sklearn.cluster import KMeans, SpectralClustering
from sklearn import preprocessing

import skfuzzy as fuzz

import f_time_frequency as tf

#%% function defination
def read_bin(fpath):
    with open(fpath, 'rb') as f:  data = f.read()
    #len(data)
    data = struct.unpack("=2901000f", data)
    data = np.array(data)
    data = data.reshape(967,3000).T
    return data


def rickerz(f=20, length=0.1, dt=0.004):
    t = np.linspace(-length/2, (length)/2, length/dt)
    y = (1-2*(np.pi**2)*(f**2)*(t**2))*np.exp(-(np.pi**2)*(f**2)*(t**2))
    return t, y


def snrz(g, f):
    '''
    To calculate the signal-to-noise_ratio.
    
    g: ground truth image
    f: noisy/restored image
    g,f: 1D or 2D array
    
    according to:
        Yangkang, Chen;
        https://en.wikipedia.org/wiki/Signal-to-noise_ratio
    '''
    g=np.array(g)
    f=np.array(f)
    
    if g.shape != f.shape: print('Dimesion of two images don''t match!')
    else:
        psnr = 20*np.log10(np.linalg.norm(g,'fro')/np.linalg.norm(g-f,'fro'))
    
    return psnr
#%
    
def seisplot_wig(s, inc=1, scale=0.8, lw=1, highlight=False, lightstep=10, 
                 figsize=(7.5,6)):
    
    nt = s.shape[0]
    xmax = s.shape[1]
    #vmax = np.max(s1)
    t = np.arange(nt)
    
    fig, ax = plt.subplots(figsize=figsize)
    for i in np.arange(0,xmax,inc):
        x1 = scale * s[:,i] + i
        
        ax.plot(x1, t, 'k', lw=lw)
        ax.fill_betweenx(t, x1, x2=i, where=x1>=i, facecolor='k', 
                          interpolate=True)
        if highlight==True:
            if i % lightstep == 0:
                ax.plot(x1, t, 'r', lw=lw*2)
                ax.fill_betweenx(t, x1, x2=i, where=x1>=i, facecolor='r', 
                              interpolate=True)
            if i==int(xmax/lightstep)*lightstep:
                ax.plot(x1, t, 'r', lw=lw*2, label='Training traces')
                ax.legend(loc='upper right' , fontsize=12)
    ax.invert_yaxis()
    ax.set_xlim(-1,xmax)
    ax.set_ylim(nt,0)
    ax.set_ylabel('Time samples', fontsize=13)
    ax.set_xlabel('Traces', fontsize=13)
    fig.tight_layout()
    
    return fig, ax

def seisplot_den(s, vmax=1):
    
    fig, ax = plt.subplots(figsize=(8,6))
    cax=ax.imshow(s,interpolation='quadric',aspect='auto',
                  vmax=vmax,vmin=-vmax)
    ax.set_ylim(s.shape[0]-1, 0)
    ax.set_xlim(0, s.shape[1]-1)
    ax.set_ylabel('Time samples', fontsize=13)
    ax.set_xlabel('Traces', fontsize=13)
    fig.colorbar(cax, shrink=0.75)
    fig.tight_layout()
    
    return fig, ax
    
    
def trace_samp(trace, tw, tinc, ws=None):
    '''
    trace sampling. v2.0
    
    trace: the trace data
    tw: time window, the window step equals to the time window
    ws: window steps
    '''
    #%
    nt = trace.shape[0]  # the number of time records of the trace
    
    
    if ws==None: ws=tw
    if nt%tw !=0: print('The rows should integer times to time window')
    
    wn = int((nt-tw)/ws + 1 ) 
    samps = np.zeros((wn,int(tw/tinc)))
    
    for i in np.arange(wn)*ws:
        win = trace[np.arange(i,i+tw,tinc)]
        samps[int(i/ws),:] = win
        
    return samps


def trace_cwt_samp(trace, dt, freq_index, tw, tinc, ws=None):
    '''
    trace sampling. v2.0
    
    trace: the trace data
    freq_index: the index of frequencies
        for example, [26,27,29,31,33,35,37,40,44,48,53,59,66,76,89]
    
    tw: time window, the window step equals to the time window
    ws: window steps
    
    for example:
        tt = trace_cwt_samp(obs, dt, freq_index, tw, dt)
    '''
    nt = trace.shape[0]  # the number of time records of the trace
    if ws==None: ws=tw
    wn = int((nt-tw)/ws + 1 ) 
    samps = np.zeros(( wn, len(freq_index), int(tw/tinc) ))
    
    cwtmatr, freqs = tf.spectrum_cwt(np.arange(nt)*dt*1000, trace, 
                            colorscale=0.5, wavelet='morl', widths=100, 
                            freqmax=100,  plot=False)
    
    for i in np.arange(wn)*ws:
        win0 = cwtmatr[freq_index,:]
                     #[20,19.4,18,17,16,15,14.2,13.2,12,11,10,9,8,7,6]
        win = win0[:,np.arange(i,i+tw,tinc)]
        samps[int(i/ws),:,:] = win
        
    return samps



def trace_samp_plot(trace, tw,tinc, nrows,ncols, 
                    ann=False, Y=[], Y_predit=[], annx=5, anny=0.15, ymax=0.2):
    '''
    plot the trace samples.
    
    trace: the trace data
    tw: time window, the window step equals to the time window
    nrows: the number of rows of the figure
    ncols: the number of colums of the figure
    
    ann: True or False,
        plot the result or not
    Y: the target result label
    Y_predict: the predicted label
    '''    
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*2,nrows*1.5), 
                            sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0.05, hspace=0.08) #left=0.08, right=0.98, 
    
    for i in np.arange(0, nrows):
        for j in np.arange(0, ncols):
            
            #if i*ncols+j > nt/tw-1: break
            
            ax = axs[i, j]
            ax.plot(np.arange(0,tw,tinc), trace[i*ncols+j,:], 'k', lw=1)
            
            if ann==True:
                ax.annotate(str(int(Y[i*ncols+j])), xy=(annx,anny), color="b")
                ax.annotate(str(int(Y_predit[i*ncols+j])), 
                            xy=(tw-annx*2,anny), color="r") 
            
            if j==0: ax.set_ylabel('Amplitude', fontsize=13)
            if i==nrows-1: ax.set_xlabel('Time samples', fontsize=13)
            ax.set_ylim(-ymax,ymax)
    #fig.tight_layout()


def trace_samp_fb_plot(isamps, tw,tinc, nrows,ncols, 
                       fb_plot=False, fb1=[],fb2=[],
                       ann=False, Y=[], Y_predit=[], 
                       annx=5, anny=0.15, ymax=0.2):
    '''
    plot the trace samples.
    
    isamps: the trace sample data
    tw: time window, the window step equals to the time window
    nrows: the number of rows of the figure
    ncols: the number of colums of the figure
    
    fb: True or False
        to plot the first-break or not
    fw: the first effective window sample
    t1,y1: the time and amplitude of the first method
    t2,y2: the time and amplitude of the second method
    
    ann: True or False,
        plot the result or not
    Y: the target result label
    Y_predict: the predicted label
    '''    
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*1.8,nrows*1.5), 
                            sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0.1, hspace=0.15, 
                        left=0.06, right=0.98, bottom=0.12, top=0.95)
    
    for i in np.arange(0, nrows):
        for j in np.arange(0, ncols):
            
            #if i*ncols+j > nt/tw-1: break
        
            ax = axs[i, j]
            ax.plot(np.arange(0,tw,tinc), isamps[i*ncols+j,:], 'k', lw=1)
            
            if fb_plot==True:
                
                fw1=int(fb1/tw)
                fw2=int(fb2/tw)
            
                t1 = fb1%tw
                t2 = fb2%tw
                
                if i*ncols+j==fw1:
                    y1=isamps[i*ncols+j,:][t1]
                    ax.scatter(t1,y1, s=120, facecolors='none', marker='^',
                               edgecolors='r', linewidths=2 )
                if i*ncols+j==fw2:
                    y2=isamps[i*ncols+j,:][t2]
                    ax.scatter(t2,y2, s=120, facecolors='none', 
                               edgecolors='b', linewidths=2 )

                if i*ncols+j==fw2:
                    ax.plot(np.arange(0,tw,tinc),isamps[i*ncols+j,:],'g',lw=2)
                if i*ncols+j==fw2+1:
                    ax.plot(np.arange(0,tw,tinc),isamps[i*ncols+j,:],'g',lw=2)
                if i*ncols+j==fw2+2:
                    ax.plot(np.arange(0,tw,tinc),isamps[i*ncols+j,:],'g',lw=2)
            
            if ann==True:
                ax.annotate(str(int(Y[i*ncols+j])), 
                            xy=(annx,anny), color="k",weight="bold")
                ax.annotate(str(int(Y_predit[i*ncols+j])), 
                            xy=(tw-annx*3,anny), color="g",weight="bold") 
            
            if j==0: ax.set_ylabel('Amplitude', fontsize=13)
            if i==nrows-1: ax.set_xlabel('Time samples', fontsize=13)
            ax.set_ylim(-ymax,ymax)
    #fig.tight_layout()
    return fig


def gather_to_label(gather, tw, tinc, thre=0.05):
    '''
    generate gather labels
    
    gather: array-like gather
    tw: time window, the window step equals to the time window
    thre: threshold to determine the labels
    '''
    
    nt = gather.shape[0]
    
    if nt%tw !=0: print('The rows should integer times to time window')
        
    labels = np.zeros((int(nt/tw), gather.shape[1]))
    labels_attr = np.zeros((int(nt/tw), gather.shape[1]))
    for i in np.arange(gather.shape[1]):
        trace_samps = trace_samp(gather[:,i], tw, tinc)
        trace_samps_attr = np.sqrt(np.mean(trace_samps**2, axis=1))
        #trace_samps_attr = np.max(np.abs(trace_samps), axis=1)
        trace_samps_class = to_classes(trace_samps_attr, threshold=thre)
        labels[:,i] = trace_samps_class
        labels_attr[:,i] = trace_samps_attr
    return labels, labels_attr


def gather_plot(gather, tw, plot_label=False, labels=[]):
    '''
    plot the gather and labels
    
    gather: array-like gather
    plot_label: True or False
    labels: labels of the samples of the gather
    '''
    
    fig, ax = plt.subplots(figsize=(9,6))
    cax=ax.imshow(gather,interpolation='quadric',aspect='auto',
                   vmax=0.2,vmin=-0.2)
    
    if plot_label==True:
        for i in np.arange(labels.shape[0]):
            for j in np.arange(labels.shape[1]): #labels.shape[1]
                if labels[i,j]==1: 
                    ax.scatter(j,i*tw+tw/2,
                                #[j,j,j], [i*tw+tw/4,i*tw+tw/2,i*tw+3*tw/4], 
                                s=25, facecolors='none', 
                                edgecolors='r', linewidths=0.2)
    
    ax.set_ylim(gather.shape[0],0)
    ax.set_xlim(0,gather.shape[1])
    ax.set_ylabel('Time samples', fontsize=13)
    ax.set_xlabel('Traces', fontsize=13)
    fig.colorbar(cax, shrink=0.75)
    fig.tight_layout()
    return fig, ax

def labels_plot(labels):
    '''
    plot the original or predicted labels of the gather.
    
    labels: labels of the samples of the gather
    savepath: path to save the pic
    '''
    cbar_binary = ListedColormap(["darkgray","yellow"]) 
    fig, ax = plt.subplots(figsize=(9,6))
    cax=ax.imshow(labels,interpolation='nearest',aspect='auto',
                   cmap=cbar_binary, vmax=1,vmin=0)
    ax.set_ylabel('Time samples', fontsize=13)
    ax.set_xlabel('Traces', fontsize=13)
    cbar = fig.colorbar(cax, shrink=0.5, ticks=[0.25,0.75])
    cbar.ax.set_yticklabels(['Non-\nWave', 'Wave'])
    fig.tight_layout()
    return fig, ax


def feature_2dplot(inter_output, samp_num, n=20):

    '''plot interlayer output 2D feature maps'''
    maximum = np.max(inter_output[samp_num,:,:,:])
    plt.figure(figsize=(10,2*n/5))
    for i in np.arange(0,n,1):
        plt.subplot(n/5,5,i+1)
        cwtmatr2 = inter_output[samp_num,i,:,:] 
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


def feature_1dplot(inter_output, samp_num, n=20):

    '''plot interlayer output 1D feature maps'''
    
    plt.figure(figsize=(12,1.75*n/5))
    for i in np.arange(0,n,1):
        plt.subplot(n/5,5,i+1)
        cwtmatr2 = inter_output[samp_num,i,0,:] 
        #maximum = 0.2       
        plt.plot(np.arange(inter_output.shape[3]), cwtmatr2, 'k')
        plt.ylim(-0.6, 0.6)
        plt.xlabel('Time samples')
        plt.ylabel('Amplitude')
#        ax = plt.gca()
#        ax.invert_yaxis()
#        ax.invert_xaxis()
#        ax.set_yticklabels([])
#        ax.set_xticklabels([])
#        ax.set_yticks([])
#        ax.set_xticks([])
    plt.tight_layout()
    #plt.colorbar()
    
        

def filter_2dplot(W, n=20):
    '''plot convolution layer 2D filters (weights)
    
    W: the weights maxtrice
    n:  the number of figures
    '''
    
    plt.figure(figsize=(10,2*n/5))
    for i in np.arange(0,n,1):
        plt.subplot(n/5,5,i+1)
        filter1 = W[:,:,0,i]       
        plt.imshow(filter1,cmap='binary')  #'RdBu'
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


def filter_1dplot(W, n=20):
    '''plot convolution layer 1D filters (weights)
    
    W: the weights maxtrice
    n:  the number of figures
    '''
    
    plt.figure(figsize=(12,1.2*n/5))
    for i in np.arange(0,n,1):
        plt.subplot(n/5,5,i+1)
        filter1 = W[0,:,0,i]
        plt.plot(filter1,'k')  #'RdBu'
        plt.xlabel('Time samples')
        plt.ylabel('Amplitude')
        #plt.ylim(10, 80)
#        ax = plt.gca()
#        ax.invert_yaxis()
#        ax.invert_xaxis()
#        ax.set_yticklabels([])
#        ax.set_xticklabels([])
#        ax.set_yticks([])
#        ax.set_xticks([])
    plt.tight_layout()
    #plt.colorbar()

def to_classes(Y_predict, threshold=0.5):
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


def pred(model, gather, g_labels, tw, tinc, 
         method, trace_step, dt, freq_index=[]):
    '''
    model: the established model, for example, MLP, CNN
    gather: the original gather
    g_labels: the labels of the gather
    tw: time window of the samples
    tinc: the time increment of samples, inside time window
    method: 'mlp', 'cnn', 'dnn', 'cwt'
    trace_step: the trace step to do prediction
    dt: the time sample interval of the gather data
    '''
    nt = gather.shape[0]
    traces = np.arange(0,gather.shape[1],trace_step)
    Y_pred = np.zeros((int(nt/tw), len(traces)))
    for i in traces:
        itrace = gather[:,i]
        
        if method=='cwt': 
            isamps = trace_cwt_samp(itrace, dt, freq_index, tw, tinc)
            isamps = isamps.reshape(-1,1,len(freq_index), int(tw/tinc))
        else:
            isamps = trace_samp(itrace, tw, tinc)
            if method=='cnn': 
                isamps = isamps.reshape(-1,1,1,int(tw/tinc))

        itrace_pred = model.predict(isamps)
        itrace_pred_class = to_classes(itrace_pred[:,1], 0.5)
        Y_pred[:,int(i/trace_step)] = itrace_pred_class
    
    '''evaluation on the result'''
    Y_pred_resh = Y_pred.reshape(Y_pred.shape[0]*Y_pred.shape[1],)
    g_labels = g_labels[:,traces]
    g_labels_resh = g_labels.reshape(Y_pred.shape[0]*Y_pred.shape[1],)
    print('The accuracy of whole gather: ',
          accuracy_score(Y_pred_resh,g_labels_resh))
    print('The F1 of whole gather: ', f1_score(Y_pred_resh,g_labels_resh))
    
    return Y_pred


def fb_pick_trace(gather, tr_num, tw, tinc, labels_pred, pred_tr_step,
                  method, n=1, c=2, w=11, nlta=20, nsta=10, a=0):
    '''
    pick the first break of one trace in the first effective wave window 
    according to the predicted labels.
    
    gather: seismic gather data
    tr_num trace number
    tw, tinc:
    
    labels_pred: the predicted labels
    pred_tr_step: the trace step of the predicted labels

    method: 'kmeans', 'fuzzy'
    
    n: the neighbouring window for clustering
    
    w: the window to calculate the attributes for clustering, odd
    c: the number of clusters
    nlta: the window for lta
    nsta: the window for sta
    a: the correction of the picked arrival time
    
    '''
#%
#    method='kmeans'
#    gather=s1
#    tr_num=14
#    labels_pred=Y_pred_cwtc
#    c=3
    
    #n = 1 # 0,1,2
    
    
    nt = gather.shape[0]
    isamps = trace_samp(gather[:, tr_num], tw, tinc)    
    
    # j: j-th sample is the first wave
    for j in np.arange(int(nt/tw)):
        if labels_pred[j, int(tr_num/pred_tr_step)]==1: break

    if j<1: n=0
    if j>=int(nt/tw): n=0
    fwtrace = isamps[np.arange(j-n,j+n+1,1),:].reshape(
                                              int(isamps.shape[1]*(2*n+1)),).T

    trace = np.array(fwtrace)
    
    '''single point attributes'''
    tr_abs = np.abs(trace).reshape(-1,1)
    
    tr_sq = trace**2
    tr_sq = tr_sq.reshape(-1,1)
    
    '''multiple point attributes'''
    adds = np.zeros(int(w/2))+(trace[0]+trace[1])/2
    adde = np.zeros(int(w/2))+(trace[-2]+trace[-1])/2
    tr = np.append(adds,trace)
    tr = np.append(tr,adde)
    
    tr_mean = np.convolve(tr, np.ones((w,))/w, mode='valid').reshape(-1,1)
    tr_absm = np.convolve(np.abs(tr), np.ones((w,))/w, 
                          mode='valid').reshape(-1,1)
    tr_sqm = np.convolve(tr**2, np.ones((w,))/w, 
                           mode='valid').reshape(-1,1)

    addlta = np.zeros(nlta-1)+(trace[0]+trace[1]+trace[2])/2
    addsta = np.zeros(nsta-1)+(trace[0]+trace[1]+trace[2])/2 
    tr_lta = np.append(addlta, trace)
    tr_sta = np.append(addsta, trace)
    
    tr_lta = np.convolve(np.abs(tr_lta), np.ones((nlta,))/nlta, 
                         mode='valid').reshape(-1,1)
    tr_sta = np.convolve(np.abs(tr_sta), np.ones((nsta,))/nsta, 
                         mode='valid').reshape(-1,1)
    tr_stalta = tr_sta/tr_lta

    '''clustering'''
    tr_comb = np.hstack((tr_sqm, tr_absm, tr_stalta)) #, tr_mean
    min_max_scaler = preprocessing.MinMaxScaler()
    tr_comb = min_max_scaler.fit_transform(tr_comb)
    
    if method=='kmeans':
        cluster = KMeans(n_clusters=c, random_state=0).fit(tr_comb)
        cntr = cluster.cluster_centers_
        clabel = cluster.labels_
    if method=='fuzzy':
        cntr, u_orig, _, _, _, _, _ = fuzz.cluster.cmeans(tr_comb.T, c=c, m=2,
                                        error=0.005,maxiter=1000,seed=201803)
        clabel = u_orig.T.argmax(axis=1)

    '''the minmum classes labels'''
    cntrindex = np.argsort(cntr[:,0], axis=0)
    cntrminindex = int(cntrindex[0])

    for z in np.arange(len(tr_abs)):
        if clabel[z] != cntrminindex: break    
    
    t_num = (j-n)*tw+z*tinc + a  # trace time sample number

#    fw = j
#    tr_t = np.arange(-tw*n,tw*(n+1),tinc).reshape(-1,1)    
#    tr_amp = trace.reshape(-1,1)
#    t = tr_t[z+a, 0]
#    y = tr_amp[z+a, 0] 
    
#    fig, ax = plt.subplots()
#    ax.plot(np.arange(-tw*n,tw*(n+1),tinc), tr, 'k')
#    ax.plot(np.arange(-tw*n,tw*(n+1)-(w-1)*tinc,tinc) +int(w/2)*tinc, tr_absm)
#    ax.plot(np.arange(-tw*n,tw*(n+1)-(w-1)*tinc,tinc) +int(w/2)*tinc, tr_sqsum)
#    ax.scatter(t, y, s=80, facecolors='none',edgecolors='b',linewidths=2)
#    ax.set_ylim(-2, 2)
#    ax.set_ylabel('Amplitude',fontsize=13)
#    ax.set_xlabel('Time sample',fontsize=13)    
#%
    return t_num



def fb_pick_wholetrace(trace, method='kmeans', plot=False, plotw=False,
                       c=2, w=11, nlta=20, nsta=10, a=0):
    '''
    trace first break picker, without window
    
    method: 'kmeans', 'fuzzy'

    c: the number of clusters
    w: the window to calculate the attributes for clustering, odd
    nlta: the window for lta
    nsta: the window for sta
    a: the correction of the picked arrival time
    '''
#    #%%
#    trace=xn
#    w = 21 # odd
#    c=2
#    method = 'fuzzy'
#    nlta = 200
#    nsta = 100
    
    trace = np.array(trace)
    
    '''single point attributes'''
    tr_abs = np.abs(trace).reshape(-1,1)
    
    tr_sq = trace**2
    tr_sq = tr_sq.reshape(-1,1)
    
    '''multiple point attributes'''
    adds = np.zeros(int(w/2))+(trace[0]+trace[1])/2
    adde = np.zeros(int(w/2))+(trace[-2]+trace[-1])/2
    tr = np.append(adds,trace)
    tr = np.append(tr,adde)
    
    tr_mean = np.convolve(tr, np.ones((w,))/w, mode='valid').reshape(-1,1)
    tr_absm = np.convolve(np.abs(tr), np.ones((w,))/w, 
                          mode='valid').reshape(-1,1)
    tr_sqm = np.convolve(tr**2, np.ones((w,))/w, 
                           mode='valid').reshape(-1,1)

    addlta = np.zeros(nlta-1)+(trace[0]+trace[1]+trace[2])/2
    addsta = np.zeros(nsta-1)+(trace[0]+trace[1]+trace[2])/2 
    tr_lta = np.append(addlta, trace)
    tr_sta = np.append(addsta, trace)
    
    tr_lta = np.convolve(np.abs(tr_lta), np.ones((nlta,))/nlta, 
                         mode='valid').reshape(-1,1)
    tr_sta = np.convolve(np.abs(tr_sta), np.ones((nsta,))/nsta, 
                         mode='valid').reshape(-1,1)
    tr_stalta = tr_sta/tr_lta

    '''clustering'''
    tr_comb = np.hstack((tr_sqm, tr_absm, tr_stalta)) #, tr_mean
    min_max_scaler = preprocessing.MinMaxScaler()
    tr_comb = min_max_scaler.fit_transform(tr_comb)
    
    if method=='kmeans':
        cluster = KMeans(n_clusters=c, random_state=0).fit(tr_comb)
        cntr = cluster.cluster_centers_
        clabel = cluster.labels_
    if method=='fuzzy':
        cntr, u_orig, _, _, _, _, _ = fuzz.cluster.cmeans(tr_comb.T, c=c, m=2,
                                        error=0.005,maxiter=1000,seed=201803)
        clabel = u_orig.T.argmax(axis=1)
    
    '''the minmum classes labels'''
    cntrindex = np.argsort(cntr[:,0], axis=0)
    cntrminindex = int(cntrindex[0])


    for z in np.arange(len(tr_abs)):
        if clabel[z] != cntrminindex:  break
    
    
    t_num = z + a  # trace time sample number
    
    if plot==True:
        fig, (ax1, ax2, ax3, ax4)= plt.subplots(4,1,figsize=(10,8))
        fig.subplots_adjust(hspace=.2)
        
        ax1.plot(np.arange(len(trace))*0.004, trace, 'k', lw=1)
        ax1.scatter(t_num*0.004, trace[t_num], s=80, facecolors='none', 
                    edgecolors='r', lw=2, marker='^')
        ax1.set_ylabel('Amplitude', fontsize=13)
        ax1.set_title('(a) Synthetic seismic', fontsize=14)     
        ax1.set_ylim(-1,2)
        if plotw==True:
            ax1.plot(np.arange(70,110)*0.004, trace[70:110], 'g', lw=2)
            nn=90 #picked arrival inside the window
            ax1.scatter(nn*0.004, trace[nn], s=100, facecolors='none', 
                    edgecolors='b', lw=2)

        ax2.plot(np.arange(len(trace))*0.004, tr_comb[:,1], 'k', lw=1)
        ax2.scatter(t_num*0.004, tr_comb[:,1][t_num], s=80, facecolors='none', 
                    edgecolors='r', lw=2, marker='^')
        ax2.set_ylabel('Amplitude', fontsize=13)
        ax2.set_title('(b) Mean absolute', fontsize=14)  
        if plotw==True:
            ax2.plot(np.arange(70,110)*0.004, tr_comb[:,1][70:110], 'g',lw=2)
            ax2.scatter(nn*0.004, tr_comb[:,1][nn], s=100, facecolors='none', 
                    edgecolors='b', lw=2)
        
        ax3.plot(np.arange(len(trace))*0.004, tr_comb[:,0], 'k', lw=1)
        ax3.scatter(t_num*0.004, tr_comb[:,0][t_num], s=80, facecolors='none', 
                    edgecolors='r', lw=2, marker='^')
        ax3.set_ylabel('Amplitude', fontsize=13)
        ax3.set_title('(c) Mean square', fontsize=14)  
        if plotw==True:
            ax3.plot(np.arange(70,110)*0.004, tr_comb[:,0][70:110], 'g',lw=2)
            ax3.scatter(nn*0.004, tr_comb[:,0][nn], s=100, facecolors='none', 
                    edgecolors='b', lw=2)
            
        ax4.plot(np.arange(len(trace))*0.004, tr_comb[:,2], 'k', lw=1)
        ax4.scatter(t_num*0.004, tr_comb[:,2][t_num], s=80, facecolors='none', 
                    edgecolors='r', lw=2, marker='^')
        ax4.set_ylabel('Amplitude', fontsize=13)
        ax4.set_xlabel('Time (s)', fontsize=13)
        ax4.set_title('(d) STA/LTA', fontsize=14)  
        if plotw==True:
            ax4.plot(np.arange(70,110)*0.004, tr_comb[:,2][70:110], 'g',lw=2)
            ax4.scatter(nn*0.004, tr_comb[:,2][nn], s=100, facecolors='none', 
                    edgecolors='b', lw=2)
        
        fig.tight_layout()
        #ax.plot(np.arange(len(trace)), tr_comb[:,0], 'y')
        #ax.plot(np.arange(len(trace)), tr_comb[:,1], 'r')
        #ax.plot(np.arange(len(trace)), tr_comb[:,2], 'b')
        #ax.plot(np.arange(len(trace)), tr_comb[:,3], 'g')
        #ax.plot(np.arange(len(tr_stalta)), tr_comb[:,3], 'g', label='mean')
        #ax1.scatter(t_num, 0, s=100, facecolors='none', edgecolors='r', lw=2)
        #ax.legend(loc='upper right' , fontsize=12)
        return fig

    #%
    else: return t_num


def fb_pick_gather(gather, tw, tinc, fb_tr_step, 
                   labels_pred, pred_tr_step, method,  
                   n=1, c=2, w=11, nlta=20, nsta=10, a=0):
    '''
    to pick the first break of the gather in the first effective wave window 
    according to the predicted labels.
    
    gather: seismic gather data
    tw, tinc:
    fb_tr_step: the trace step to pick first break
    
    labels_pred: the predicted labels
    pred_tr_step: the trace step of the predicted labels
    
    method: 'kmeans', 'fuzzy'
    
    w: the window to calculate the attributes for clustering, odd
    c: the number of clusters
    '''
    fb = np.zeros((int((gather.shape[1]-1)/fb_tr_step)+1,2))
    
    for i in np.arange(0, gather.shape[1], fb_tr_step):
#            ifb, _, _, _, _= fb_pick_trace_fz(gather, i, tw, tinc, fb_tr_step,
#                                      labels_pred, pred_tr_step, w=w, c=c)
        ifb = fb_pick_trace(gather, i, tw, tinc, 
                                     labels_pred, pred_tr_step, method=method,  
                                     n=n, c=c, w=w, nlta=nlta, nsta=nsta, a=a)
        fb[int(i/fb_tr_step),:] = [i,ifb]
    
    return fb


def fb_pick_gather_wholetrace(gather, fb_tr_step, method,
                              c=2, w=11, nlta=20, nsta=10, a=0):
    '''
    to pick the first break of the gather in the whole trace.
    
    gather: seismic gather data
    
    method: 'kmeans', 'fuzzy'
    
    fb_tr_step: the trace step to pick first break

    w: the window to calculate the attributes for clustering, odd
    c: the number of clusters
    '''
    fb = np.zeros((int((gather.shape[1]-1)/fb_tr_step)+1,2))
    
    for i in np.arange(0, gather.shape[1], fb_tr_step):

        tr = gather[:,i]
        ifb = fb_pick_wholetrace(tr, method, 
                                 w=w, c=c, nlta=nlta, nsta=nsta, a=a)

        fb[int(i/fb_tr_step),:] = [i,ifb]
    
    return fb







