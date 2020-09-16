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

#%% Package import

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier 
from sklearn.metrics import accuracy_score, f1_score
from sklearn import preprocessing

from keras.utils import np_utils
from keras.optimizers import Adam
from keras.models import Model
from keras.models import load_model

import f_nnclassifier as classifier
import f_time_frequency as tf
import f



#%% multi-channel synthetic data
''' multi-channel synthetic data (clean) '''
s0 = np.array(pd.read_csv("data/syn_data.csv", header=None))
#fig, ax = plt.subplots(figsize=(10,8))
#ax.imshow(s1,interpolation='quadric', aspect='auto',vmax=0.2,vmin=-0.2)

s0 = s0[350:350+360, :][:, 40:101]
min_max_scaler = preprocessing.MaxAbsScaler()
s0 = min_max_scaler.fit_transform(s0)
#fig, ax = plt.subplots(figsize=(10,8))
#ax.imshow(s1,interpolation='quadric', aspect='auto')

s0fb = f.fb_pick_gather_wholetrace(s0,fb_tr_step=1,method='kmeans',w=21, c=2)

fig, ax = f.seisplot_wig(s0)
#ax.scatter(s0fb[:,0],s0fb[:,1], s=50,facecolors='none',edgecolors='b',lw=1)
#fig.savefig('fig/s1.pdf', dpi=200)


#%
''' multi-channel synthetic data (noisy) '''
np.random.seed(20180320)  #20180320
s0n = s0 + np.random.randn(s0.shape[0],s0.shape[1])*0.2

f.fb_pick_wholetrace(s0n[:,14], method='kmeans', plot=True,
                       c=2, w=21, nlta=20, nsta=10, a=0)

#%
ssnr = f.snrz(s0,s0n)
print('The signal to noise ratio is:', ssnr)

dt = 0.001
t = np.arange(360)*dt*1000
cwtmatr, freqs =tf.spectrum_cwt(t, s0n[:,13], wavelet='morl', 
                 widths=100, cmap='RdBu', colorscale=1, contour_dvision=41, 
                 freqmax=100, figsize=(12,3), plot=True, ymin=6, ymax=90)

s0nfb = f.fb_pick_gather_wholetrace(s0n,fb_tr_step=1,method='kmeans',
                                    c=2, w=21, nlta=20, nsta=10, a=0)

fig, ax = f.seisplot_wig(s0n, lw=0.8)
ax.scatter(s0nfb[:,0],s0nfb[:,1], s=50,facecolors='none',edgecolors='r',lw=2)
ax.set_ylabel('Time (ms)', fontsize=13)
#fig.savefig('fig/s1n_snr=%.2f.png'%ssnr, dpi=200)

#%% parameter defination and label generate

dt = 0.001  # time interval
freq_index = [4,5,6,7,8,9,10,11,13,16,19,24,32,50,98]
freq = freqs[freq_index,]

tw = 20   # time window
tinc = 1  # time increment of samples, inside time window
nt = s0.shape[0]  # the number of time records

ncols = 6
nrows = int(nt/tw/ncols)

testsize = 0.3
seed = 20180319
channel = 1 
length = int(tw/tinc)
classes = 2

cbar_binary = ListedColormap(["darkgray","yellow"]) 



s1_labels, labels_attr = f.gather_to_label(s0, tw, tinc, thre=0.35)

fig,ax = f.labels_plot(s1_labels)
ax.set_yticks(np.arange(0,nt/tw,2.5)-0.5)
ax.set_yticklabels( ((np.arange(0,nt/tw,2.5)-0.5)*tw+tw/2).astype(int) )
ax.set_ylabel('Time (ms)', fontsize=13)
#fig.savefig('fig/s0n_labels.pdf'%ssnr, dpi=200)

s1 = s0n


#%% training dataset for MLP, CNN, and CWT-CNN

traces_step = 4
traces_train = np.arange(0, s1.shape[1], traces_step)

X = np.zeros((int(len(traces_train)*nt/tw), int(tw/tinc) ))
Xcwt = np.zeros((int(len(traces_train)*nt/tw), len(freq_index), int(tw/tinc)))
Y = np.zeros((int(len(traces_train)*nt/tw), ))
for i in traces_train:
    itrace = s1[:,i]

    itrace_samps = f.trace_samp(itrace, tw, tinc)
    X[int(i/traces_step*nt/tw):int((i/traces_step+1)*nt/tw),:] = itrace_samps
    
    isamps_cwt = f.trace_cwt_samp(itrace, dt, freq_index, tw, tinc)
    Xcwt[int(i/traces_step*nt/tw):int((i/traces_step+1)*nt/tw),:,:]=isamps_cwt
    
    ilabel = s1_labels[:,i]
    Y[int(i/traces_step*nt/tw):int((i/traces_step+1)*nt/tw),] = ilabel

'''plot the selected traces'''
fig,ax = f.seisplot_wig(s1, scale=0.8, lw=0.5, 
                        highlight=True, lightstep=traces_step)
#ax.set_yticks(np.arange(0,nt,50))
#ax.set_yticklabels(np.arange(0,nt,50)*dt )
ax.set_ylabel('Time (ms)', fontsize=13)
#fig.savefig('fig/s0n_partition_snr=%.2f.pdf'%ssnr, dpi=200)

''' data partitioning'''
X_train, X_test, Y_train, Y_test = \
             train_test_split(X, Y, test_size=testsize, random_state=seed)
Xcwt_train, Xcwt_test, Y_train, Y_test = \
             train_test_split(Xcwt, Y, test_size=testsize, random_state=seed)
Y_train = np_utils.to_categorical(Y_train, num_classes=classes)
Y_test = np_utils.to_categorical(Y_test, num_classes=classes) 

Xcwt_train = Xcwt_train.reshape(-1, channel, len(freq_index), length)
Xcwt_test = Xcwt_test.reshape(-1, channel, len(freq_index), length)  

#%% CWT-CNN

modelcwtc, histcwtc = classifier.cwtcnn(Xcwt_train,Xcwt_test,Y_train,Y_test, 
                   activation='relu', neurons=50, hflayer=2, drop=0.3, 
                   optimizer=Adam(lr=0.0001),  #'SGD' Adam(lr=0.0001)
                   batchsize=40, epoches=1000,
                   reduce_threshold=0.005, lr_patience=10, lr_reduce=0.2,
                   stop_threshold=0.005, stop_patience=20,
                   filters=40, kernel_size=3, strides=2)
classifier.plot(histcwtc)

Y_train_cwtc = modelcwtc.predict(Xcwt_train)
Y_test_cwtc = modelcwtc.predict(Xcwt_test)

Y_train_cwtc_class = f.to_classes(Y_train_cwtc[:,1], 0.5)
Y_test_cwtc_class = f.to_classes(Y_test_cwtc[:,1], 0.5)

print('The accuracy of training set: ', 
      accuracy_score(Y_train[:,1],Y_train_cwtc_class))
print('The accuracy of testing set:', 
      accuracy_score(Y_test[:,1], Y_test_cwtc_class))
print('The F1 of training set: ', f1_score(Y_train[:,1],Y_train_cwtc_class))
print('The F1 of testing set:', f1_score(Y_test[:,1],Y_test_cwtc_class))


#%% CWT-CNN  --  analysis

layer_list = list([(layer.name, layer) for layer in modelcwtc.layers])
#modelcwtc.summary()

''' calculate the interlayer output of CWT-CNNs'''
layer_name = layer_list[0][0]  #'block5_conv3' 'conv2d_1'
inter_layer_model = Model(inputs=modelcwtc.input,
                          outputs=modelcwtc.get_layer(layer_name).output)

tr_num = 13
obs_samps_cwt = f.trace_cwt_samp(s1[:,tr_num], dt, freq_index, tw, tinc)
obs_samps_cwt = obs_samps_cwt.reshape(-1,1,len(freq_index),length)
inter_output = inter_layer_model.predict(obs_samps_cwt) 

#%
'''plot the input maps '''
samp_num = 6

cwtmatr = obs_samps_cwt[samp_num, 0,:,:]
cwtmatrT = cwtmatr.T
maximum = np.max(cwtmatrT)
cwtmatrT = np.clip(cwtmatrT, -maximum, maximum, cwtmatrT)
levels = np.linspace(-maximum, maximum, 41)
fig, ax = plt.subplots(figsize=(5,4))
cax=ax.contourf(freq,np.arange(0,tw,tinc), cwtmatrT,levels=levels,cmap='RdBu')
#plt.ylim(10, 80)
ax.set_xlim(6, 40)
ax.invert_yaxis()
ax.set_xlabel('Frequency(Hz)')
ax.set_ylabel('Time samples')
fig.colorbar(cax)

'''plot interlayer output feature maps'''
f.feature_2dplot(inter_output, samp_num, n=20)

'''plot convolution layer filters (weights)'''
wb = modelcwtc.layers[0].get_weights()
w = wb[0]
f.filter_2dplot(w, n=20)


#%% CNN

#import f_nnclassifier as classifier

Xcnn_train = X_train.reshape(-1, channel, 1, length)
Xcnn_test = X_test.reshape(-1, channel, 1, length) 

modelcnn, histcnn = classifier.cnn(Xcnn_train, Xcnn_test, Y_train, Y_test,
                                   activation='relu', neurons=100, hflayer=2,
                                   drop=0.3, optimizer=Adam(lr=0.0001), 
                                   batchsize=80, epoches=1000,
                                   filters=40, kernel_size=(1,5), strides=2)
classifier.plot(histcnn)
#modelcnn.save('out/model_cnn_n03_201804.h5')

Y_train_cnn = modelcnn.predict(Xcnn_train)
Y_test_cnn = modelcnn.predict(Xcnn_test)

Y_train_cnn_class = f.to_classes(Y_train_cnn[:,1], 0.5)
Y_test_cnn_class = f.to_classes(Y_test_cnn[:,1], 0.5)

print('The accuracy of training set: ', 
      accuracy_score(Y_train[:,1],Y_train_cnn_class))
print('The accuracy of testing set:', 
      accuracy_score(Y_test[:,1], Y_test_cnn_class))
print('The F1 of training set: ', f1_score(Y_train[:,1],Y_train_cnn_class))
print('The F1 of testing set:', f1_score(Y_test[:,1],Y_test_cnn_class))

#%% CNN  --  analysis
# modelcnn = load_model('out/model_cnn_n06_201804.h5')
# modelcnn.summary()

#layer_list = list([(layer.name, layer) for layer in modelcnn.layers])

tr_num = 13
samp_num = 6
 
'''plot the input maps (noisy) '''
samps = f.trace_samp(s1[:,tr_num], tw, tinc)
samps = samps.reshape(-1, channel, 1, length)
matr = samps[samp_num, 0, 0,:]

fig, ax = plt.subplots(figsize=(5.5,4))
cax=ax.plot(np.arange(0,tw,tinc), matr, 'k')
ax.set_ylim(-1.3, 1.3)
ax.set_xlabel('Time samples', fontsize=13)
ax.set_ylabel('Amplitude', fontsize=13)
ax.set_xticks([0,5,10,15,20])
ax.set_xticklabels([0,5,10,15,20])
ax.set_yticks([-1,-0.5,0,0.5,1])
ax.set_yticklabels([-1,-0.5,0,0.5,1])
#fig.savefig('fig/sampn', dpi=200)

'''plot the input maps (clean) '''
samps0 = f.trace_samp(s0[:,tr_num], tw, tinc)
samps0 = samps0.reshape(-1, channel, 1, length)
matr0 = samps0[samp_num, 0, 0,:]

fig, ax0 = plt.subplots(figsize=(5.5,4))
ax0.plot(np.arange(0,tw,tinc), matr0, 'k')
ax0.set_ylim(-1.3, 1.3)
ax0.set_xlabel('Time samples', fontsize=13)
ax0.set_ylabel('Amplitude', fontsize=13)
ax0.set_xticks([0,5,10,15,20])
ax0.set_xticklabels([0,5,10,15,20])
ax0.set_yticks([-1,-0.5,0,0.5,1])
ax0.set_yticklabels([-1,-0.5,0,0.5,1])
#fig.savefig('fig/samp', dpi=200)


'''plot the convolution layer output maps '''
conv_name = modelcnn.layers[0].name
conv_model = Model(inputs=modelcnn.input,
                          outputs=modelcnn.get_layer(conv_name).output)
conv_output = conv_model.predict(samps)

f.feature_1dplot(conv_output, samp_num, n=20)
#plt.savefig('fig/conv_feature', dpi=200)


'''plot the pooling layer output maps '''
pool_name = modelcnn.layers[1].name
pool_model = Model(inputs=modelcnn.input,
                   outputs=modelcnn.get_layer(pool_name).output)
pool_output = pool_model.predict(samps)

f.feature_1dplot(pool_output, samp_num, n=20)
#plt.savefig('fig/pool_feature', dpi=200)


'''plot convolution filters (weights)'''
wb1 = modelcnn.layers[0].get_weights()
w1 = wb1[0]
f.filter_1dplot(w1, n=20)
#plt.savefig('fig/conv_filter', dpi=200)

#%% DNN
"""
modeldnn, histdnn = classifier.dnn(X_train, X_test,Y_train, Y_test,
                       activation='relu', neurons=100, drop=0.2, hflayer=3,
                       optimizer=Adam(lr=0.0001), batchsize=40, epoches=1000)
classifier.plot(histdnn)
# model1d.save('model/model_cnn1d_20180209.h5')

Y_train_dnn = modeldnn.predict(X_train)
Y_test_dnn = modeldnn.predict(X_test)

Y_train_dnn_class = f.to_classes(Y_train_dnn[:,1], 0.5)
Y_test_dnn_class = f.to_classes(Y_test_dnn[:,1], 0.5)

print('The accuracy of training set: ', 
      accuracy_score(Y_train[:,1],Y_train_dnn_class))
print('The accuracy of testing set:', 
      accuracy_score(Y_test[:,1], Y_test_dnn_class))
print('The F1 of training set: ', f1_score(Y_train[:,1],Y_train_dnn_class))
print('The F1 of testing set:', f1_score(Y_test[:,1],Y_test_dnn_class))
"""
#%% MLP

mlp = MLPClassifier(random_state=1115, verbose=0, 
                   hidden_layer_sizes=(100,), alpha=0.01, max_iter=1000) 
mlp.fit(X_train, Y_train)
Y_train_mlp = mlp.predict(X_train)
Y_test_mlp = mlp.predict(X_test)

Y_train_mlp_class = f.to_classes(Y_train_mlp[:,1], 0.5)
Y_test_mlp_class = f.to_classes(Y_test_mlp[:,1], 0.5)

print('The accuracy of training set: ', 
      accuracy_score(Y_train[:,1], Y_train_mlp_class))
print('The accuracy of testing set:', 
      accuracy_score(Y_test[:,1], Y_test_mlp_class))
print('The F1 of training set: ', f1_score(Y_train[:,1],Y_train_mlp_class))
print('The F1 of testing set:', f1_score(Y_test[:,1],Y_test_mlp_class))

#%% Classification results on whole line

pred_tr_step = 1

'''prediction on gather'''
print('CWT-CNN:')
Y_pred_cwtc = f.pred(modelcwtc, s1,s1_labels, tw,tinc,method='cwt',
                     trace_step=pred_tr_step, dt=dt, freq_index=freq_index)
print('CNN:')
Y_pred_cnn = f.pred(modelcnn, s1,s1_labels, tw,tinc, method='cnn',
                    trace_step=pred_tr_step, dt=dt )
#print('DNN:')
#Y_pred_dnn = f.pred(modeldnn, s1,s1_labels, tw,tinc, method='dnn',
#                    trace_step=pred_tr_step, dt=dt)
print('MLP: ')
Y_pred = f.pred(mlp, s1,s1_labels, tw,tinc, method='mlp',
                trace_step=pred_tr_step, dt=dt )

#%
'''plot the result on gather'''
fig, ax = f.labels_plot(Y_pred_cwtc)
ax.set_yticks(np.arange(0,int(nt/tw),2.5)-0.5)
ax.set_yticklabels( ( (np.arange(0,int(nt/tw),2.5)-0.5)*tw+tw/2 )*dt )
ax.set_ylabel('Time (s)', fontsize=13)
#fig.savefig('fig/s0n_labels_cwtcnn_snr=%.2f.pdf'%ssnr, dpi=200)

fig, ax = f.labels_plot(Y_pred_cnn)
ax.set_yticks(np.arange(0,nt/tw,2.5)-0.5)
ax.set_yticklabels( ((np.arange(0,nt/tw,2.5)-0.5)*tw+tw/2).astype(int) )
ax.set_ylabel('Time (ms)', fontsize=13)
#fig.savefig('fig/s0n_labels_cnn_snr=%.2f.pdf'%ssnr, dpi=200)

#fig, ax = f.labels_plot(Y_pred_dnn)
#ax.set_yticks(np.arange(0,int(nt/tw),2.5)-0.5)
#ax.set_yticklabels( ( (np.arange(0,int(nt/tw),2.5)-0.5)*tw+tw/2 )*dt )
#ax.set_ylabel('Time (s)', fontsize=13)
#fig.savefig('fig/s0n_labels_dnn_snr=%.2f.pdf'%ssnr, dpi=200)

fig, ax = f.labels_plot(Y_pred)
ax.set_yticks(np.arange(0,nt/tw,2.5)-0.5)
ax.set_yticklabels( ((np.arange(0,nt/tw,2.5)-0.5)*tw+tw/2).astype(int) )
ax.set_ylabel('Time (ms)', fontsize=13)
#fig.savefig('fig/s0n_labels_mlp_snr=%.2f.pdf'%ssnr, dpi=200)
#%% Classification results on one trace
import f

tr_num = 31

s1_samps = f.trace_samp(s1[:,tr_num], tw, tinc)
f.trace_samp_fb_plot(s1_samps, tw, tinc, nrows, ncols, ann=True, 
                  Y=s1_labels[:,tr_num], 
                  Y_predit=Y_pred_cnn[:,int(tr_num/pred_tr_step)],
                  annx=1, anny=1.5, ymax=2.5)

#%% First-break pick -- gather   ##############
import time
t0 = time.time()


fb_tr_step = 1
w = 21  # window for attribute calculation
c = 2  # clusters
n = 1  # the neighbouring window for clustering
a = 0  # correction for picked arrival time

t0 = time.time()
'''1-wave classification, 2-clustering inside the first effective wave'''
fb1 = f.fb_pick_gather(s1, tw, tinc, fb_tr_step, Y_pred_cnn, 
                       pred_tr_step, method='kmeans', w=w, c=c, n=n, a=a)  
t1 = time.time()
print('time to pick arrivals using CNN-KC:', t1-t0)
#fb2 = f.fb_pick_gather(s1, tw, tinc, fb_tr_step, Y_pred_cnn, 
#                       pred_tr_step, method='fuzzy', w=w, c=c, n=n) 

'''1-clustering with one whole trace data'''
fb3 = f.fb_pick_gather_wholetrace(s1, fb_tr_step, method='kmeans', w=w, c=c) 
#fb4 = f.fb_pick_gather_wholetrace(s1, fb_tr_step, method='fuzzy', w=w, c=c) 
t2 = time.time()
print('time to pick arrivals using KC:', t2-t1)

fig, ax = f.seisplot_wig(s1, lw=0.3)
ax.scatter(fb1[:,0],fb1[:,1], s=80, 
            facecolors='none',edgecolors='b', lw=2, label='CNN-KC')
#ax.scatter(fb2[:,0],fb2[:,1], s=50, marker='D',
#            facecolors='none',edgecolors='b', lw=1, label='fuzzy')
ax.scatter(fb3[:,0],fb3[:,1], s=60, marker='^',
            facecolors='none',edgecolors='r', lw=2, label='KC')
#ax.scatter(fb4[:,0],fb4[:,1], s=60, marker='s',
#            facecolors='none',edgecolors='r', lw=1, label='fuzzy-wholetrace')
ax.legend(loc='lower right' , fontsize=12)
#ax.set_yticks(np.arange(0,nt,50))
#ax.set_yticklabels(np.arange(0,nt,50)*dt )
ax.set_ylabel('Time (ms)', fontsize=13)
#fig.savefig('fig/s0n_fb_snr=%.2f.pdf'%ssnr, dpi=200)

err_cnnkmeans = np.sum(np.abs(fb1[:,1]-s0fb[:,1]))
err_kmeans = np.sum(np.abs(fb3[:,1]-s0fb[:,1]))
print('The error of CNN-kmeans: ', err_cnnkmeans)
print('The error of Kmeans:', err_kmeans)


#%% First-break pick -- one trace

i=21  #13
isamps = f.trace_samp(s1[:,i], tw, tinc)

'''1-wave classification, 2-clustering inside the first effective wave'''
ifb1 = f.fb_pick_trace(s1, i, tw, tinc, Y_pred_cnn, pred_tr_step,
                                    method='kmeans', w=w, c=c, n=n, a=a)
#ifb2, fw2, t2, y2 = f.fb_pick_trace(s1, i, tw, tinc, Y_pred_cnn, pred_tr_step,
#                                    method='fuzzy', w=w, c=c, n=n)

#'''1-wave classification, 2-clustering with all the first effective waves
#with fuzzy clustering'''
#ifb2, fw2, t2, y2, cntr = f.fb_pick_trace_allfw(s1, i, tw, tinc, fb_tr_step,
#                                          Y_pred_cwtc, pred_tr_step, w=w, c=c)

'''1-clustering using one whole trace data with fuzzy clustering'''
ifb3 = f.fb_pick_wholetrace(s1[:,i], method='kmeans', w=w, c=c)
#ifb4 = f.fb_pick_wholetrace(s1[:,i], method='fuzzy', w=31, c=c)

#%% 
import f
fig = f.trace_samp_fb_plot(isamps, tw,tinc, nrows=3, ncols=6, fb_plot=True, 
                     fb1=ifb3, fb2=ifb1, ann=True, Y=s1_labels[:,i], 
                     Y_predit=Y_pred_cnn[:,int(i/pred_tr_step)],
                     annx=1, anny=2, ymax=3)
#fig.savefig('fig/s0n_samps_fb_n00.pdf', dpi=200)

#%% 
#fig, ax = plt.subplots()
#ax.plot(np.arange(0,tw,tinc),isamps[fw1,:], 'k')
#ax.scatter(t1,y1, s=50, facecolors='none',edgecolors='r',linewidths=2)
#ax.scatter(t2,y2, s=50, facecolors='none',edgecolors='b',linewidths=2)
#ax.scatter(t3,0 , s=50, facecolors='none',edgecolors='y',linewidths=2)
#ax.set_ylim(-0.2,0.2)
#ax.set_ylabel('Amplitude',fontsize=13)
#ax.set_xlabel('Time sample',fontsize=13)

fig, ax = plt.subplots(figsize=(13,3))
ax.plot(np.arange(nt), s1[:,i], 'k',lw=0.8)
ax.scatter(ifb1,0, s=80, facecolors='none',edgecolors='b',lw=2)
#ax.scatter(ifb2,0, s=50, marker='D', facecolors='none',edgecolors='b',lw=2)
ax.scatter(ifb3,0, s=80, marker='^', facecolors='none',edgecolors='r',lw=2)
#ax.scatter(ifb4,0, s=50, marker='s', facecolors='none',edgecolors='y',lw=2)
ax.set_ylabel('Amplitude', fontsize=13)
ax.set_xlabel('Time samples', fontsize=13)
fig.tight_layout()
ax.legend(loc='upper right' , fontsize=12)
#fig.savefig('fig/s1n_trace_fb', dpi=200)









