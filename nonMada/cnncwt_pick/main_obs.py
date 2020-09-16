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

import f_nnclassifier as classifier
import f_time_frequency as tf
import f

#%% Data import and parameter defination

'''import observed and synthetic data'''
#Syn = f.read_bin(fpath='data/syn2d_3000_967_fhi10.bin')
Obso = f.read_bin(fpath='data/obs2d_3000_967_fhi10.bin')
#Obs = pd.read_csv("data/micro.csv", header=None)

# add normally random noise
np.random.seed(20180319)
Obs = Obso  + np.random.randn(Obso.shape[0],Obso.shape[1])*0.04

ssnr = f.snrz(Obso,Obs)
print('The signal to noise ratio is:', ssnr)

nt = Obs.shape[0]  # the number of time records of a trace
dt = 0.0015  # time interval

t = np.arange(nt)*dt*1000
cwtmatr, freqs =tf.spectrum_cwt(t, Obs[:,260], wavelet='morl', #'mexh''morl'
                 widths=100, cmap='RdBu', 
                 colorscale=1, contour_dvision=41, freqmax=200,
                 figsize=(12,3), plot=True,ymin=5, ymax=30)

tw = 50  # time window size
ts = tw  # time window steps
freq_index = [26,27,29,31,33,35,37,40,44,48,53,59,66,76,89]
freq = freqs[freq_index,]

ncols = 6
nrows = int(nt/tw/ncols)

testsize = 0.2
seed = 20180319
channel = 1
tinc = 2  # time increment of samples, inside time window
length = int(tw/tinc)
classes = 2

cbar_binary = ListedColormap(["darkgray","yellow"]) 

tr_num = 260  # trace number

#%%
Obsfb = f.fb_pick_gather_wholetrace(Obs,fb_tr_step=10,method='kmeans',w=21,c=2)
fig, ax = f.seisplot_wig(Obs, inc=10, scale=30, lw=0.2, )
#ax.scatter(Obsfb[:,0],Obsfb[:,1], s=30,facecolors='none',edgecolors='r',lw=1)
ax.set_yticks(np.arange(0,nt+1,400))
ax.set_yticklabels(np.arange(0,nt+1,400)*dt)
ax.set_xticks(np.arange(0,967,100))
ax.set_xticklabels(np.arange(0,967,100)/10)
#fig.savefig('fig/obs', dpi=200)

#%% samples generate for one trace

#obs = Obs[:,tr_num]
#fig, axe = plt.subplots(figsize=(12,3))
#axe.plot(np.arange(len(obs)), obs, 'k', lw=1)
#axe.set_ylabel('Amplitude', fontsize=13)
#axe.set_xlabel('Time samples', fontsize=13)

'''generate samples for one trace'''
#obs_samps = f.trace_samp(obs, tw, tinc)
#obs_samps_cwt = f.trace_cwt_samp(obs, dt, tw, tinc)
#
#f.trace_samp_plot(obs_samps, tw, tinc, nrows, ncols)

#%% plot the gather and labels

'''generate labels for gather'''
#f.gather_plot(Obs, tw)

#Obs_labels, labels_attr_rms = f.gather_to_label(Obs, tw, tinc, thre=0.035)

#Obs_labels = labels
#np.savetxt("out/Obs_labels_win50.csv", Obs_labels, delimiter=',')
Obs_labels = pd.read_csv("out/Obs_labels_win50_m.csv", header=None)
Obs_labels = np.array(Obs_labels)

#f.gather_plot(Obs, tw, savepath='fig/gather_labels.png', 
#            plot_label=True, labels=Obs_labels)

fig,ax = f.labels_plot(Obs_labels)
ax.set_yticks(np.arange(0,nt/tw+1,400/tw))
ax.set_yticklabels(np.arange(0,nt+1,400)*dt)
ax.set_xticks(np.arange(0,967,100))
ax.set_xticklabels((np.arange(0,967,100)/10).astype(int))
#fig.savefig('fig/obs_label', dpi=200)
#%% training dataset for MLP, CNN, and CWT-CNN
traces_step = 80
traces_train = np.arange(0, 967, traces_step)

X = np.zeros((int(len(traces_train)*nt/tw), int(tw/tinc) ))
Xcwt = np.zeros((int(len(traces_train)*nt/tw), 15, int(tw/tinc) ))
Y = np.zeros((int(len(traces_train)*nt/tw), ))
for i in traces_train:
    itrace = Obs[:,i]

    itrace_samps = f.trace_samp(itrace, tw, tinc)
    X[int(i/traces_step*nt/tw):int((i/traces_step+1)*nt/tw),:] = itrace_samps
    
    isamps_cwt = f.trace_cwt_samp(itrace, dt, freq_index, tw, tinc)
    Xcwt[int(i/traces_step*nt/tw):int((i/traces_step+1)*nt/tw),:,:]=isamps_cwt
    
    ilabel = Obs_labels[:,i]
    Y[int(i/traces_step*nt/tw):int((i/traces_step+1)*nt/tw),] = ilabel

'''plot the selected traces'''
fig, ax = f.seisplot_wig(Obs, inc=10, scale=30, lw=0.2, 
                         highlight=True, lightstep=traces_step)
ax.set_yticks(np.arange(0,nt+1,400))
ax.set_yticklabels(np.arange(0,nt+1,400)*dt)
ax.set_xticks(np.arange(0,967,100))
ax.set_xticklabels((np.arange(0,967,100)/10).astype(int))
#fig.savefig('fig/obs_partition', dpi=200)

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
"""
modelcwtc, histcwtc = classifier.cwtcnn(Xcwt_train,Xcwt_test,Y_train,Y_test, 
                   activation='relu', neurons=50, hflayer=2, drop=0.3, 
                   optimizer=Adam(lr=0.0001),  #'SGD' Adam(lr=0.0001)
                   batchsize=40, epoches=1000,
                   reduce_threshold=0.005, lr_patience=10, lr_reduce=0.2,
                   stop_threshold=0.005, stop_patience=20,
                   filters=40, kernel_size=5, strides=2)
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

tr_num = 260
obs_samps_cwt = f.trace_cwt_samp(Obs[:,tr_num], dt, freq_index, tw, tinc)
obs_samps_cwt = obs_samps_cwt.reshape(-1,1,len(freq_index),length)
inter_output = inter_layer_model.predict(obs_samps_cwt) 

#%
'''plot the input maps '''
samp_num = 28

cwtmatr = obs_samps_cwt[samp_num, 0,:,:]
cwtmatrT = cwtmatr.T
maximum = np.max(cwtmatrT)
cwtmatrT = np.clip(cwtmatrT, -maximum, maximum, cwtmatrT)
levels = np.linspace(-maximum, maximum, 41)
fig, ax = plt.subplots(figsize=(5,4))
cax=ax.contourf(freq,np.arange(0,tw,tinc), cwtmatrT,levels=levels,cmap='RdBu')
#plt.ylim(10, 80)
ax.set_xlim(6, 20)
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



#%% CWT-DNN

modelcwtd, histcwtd = classifier.cwtdnn(Xcwt_train,Xcwt_test,Y_train,Y_test, 
                   activation='relu', neurons=50, hflayer=3, drop=0.5, 
                   optimizer=Adam(lr=0.0001),  #'SGD' Adam(lr=0.0001)
                   batchsize=40, epoches=1000,
                   reduce_threshold=0.005, lr_patience=10, lr_reduce=0.2,
                   stop_threshold=0.005, stop_patience=20)
classifier.plot(histcwtd)

Y_train_cwtd = modelcwtd.predict(Xcwt_train)
Y_test_cwtd = modelcwtd.predict(Xcwt_test)

Y_train_cwtd_class = f.to_classes(Y_train_cwtd[:,1], 0.5)
Y_test_cwtd_class = f.to_classes(Y_test_cwtd[:,1], 0.5)

print('The accuracy of training set: ', 
      accuracy_score(Y_train[:,1],Y_train_cwtd_class))
print('The accuracy of testing set:', 
      accuracy_score(Y_test[:,1], Y_test_cwtd_class))
print('The F1 of training set: ', f1_score(Y_train[:,1],Y_train_cwtd_class))
print('The F1 of testing set:', f1_score(Y_test[:,1],Y_test_cwtd_class))

"""
#%% CNN

Xcnn_train = X_train.reshape(-1, channel, 1, length)
Xcnn_test = X_test.reshape(-1, channel, 1, length) 

modelcnn, histcnn = classifier.cnn(Xcnn_train, Xcnn_test, Y_train, Y_test,
                                   activation='relu', neurons=100, hflayer=2,
                                   drop=0.3, optimizer=Adam(lr=0.0001), 
                                   batchsize=80, epoches=1000,
                                   filters=40, kernel_size=(1,5), strides=2, 
                                   poolstride=2)
classifier.plot(histcnn)
#model1d.save('model/model_cnn1d_20180209.h5')

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
#modelcnn = load_model('out/model_cnn_micro_201804.h5')
# modelcnn.summary()

#layer_list = list([(layer.name, layer) for layer in modelcnn.layers])
tr_num = 260
samp_num = 28
 
'''plot the input maps (noisy) '''
samps = f.trace_samp(Obs[:,tr_num], tw, tinc)
samps = samps.reshape(-1, channel, 1, length)
matr = samps[samp_num, 0, 0,:]

fig, ax = plt.subplots(figsize=(5.5,4))
cax=ax.plot(np.arange(0,tw,tinc), matr, 'k')
#ax.set_ylim(-0.5, 0.5)
ax.set_xlabel('Time samples', fontsize=13)
ax.set_ylabel('Amplitude', fontsize=13)
#ax.set_xticks([0,5,10])
#ax.set_xticklabels([0,5,10])
ax.set_yticks([-0.5,0,0.5])
#fig.savefig('fig/sampn', dpi=200)


'''plot the convolution layer output maps '''
conv_name = modelcnn.layers[0].name
conv_model = Model(inputs=modelcnn.input,
                          outputs=modelcnn.get_layer(conv_name).output)
conv_output = conv_model.predict(samps)

f.feature_1dplot(conv_output, samp_num, n=20)
#plt.savefig('fig/conv_feature', dpi=200)

'''plot the convolution layer output maps '''
conv_name = modelcnn.layers[1].name
conv_model = Model(inputs=modelcnn.input,
                          outputs=modelcnn.get_layer(conv_name).output)
conv_output = conv_model.predict(samps)

f.feature_1dplot(conv_output, samp_num, n=20)
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


#%% Results on whole line

#fig, axe = plt.subplots()
#axe.plot(np.arange(len(Y_test[:,1])), Y_test[:,1], 'ro')
#axe.plot(np.arange(len(Y_test[:,1])), Y_test_mlp_class, 'b*')
#axe.set_xlim(200,250)

pred_tr_step = 10

'''prediction on gather'''
print('MLP: ')
Y_pred = f.pred(mlp, Obs,Obs_labels, tw,tinc, method='mlp',
              trace_step=pred_tr_step, dt=dt )
#print('DNN:')
#Y_pred_dnn = f.pred(modeldnn, Obs,Obs_labels, tw,tinc, method='dnn',
#                  trace_step=pred_tr_step, dt=dt )
print('CNN:')
Y_pred_cnn = f.pred(modelcnn, Obs,Obs_labels, tw,tinc, method='cnn',
                  trace_step=pred_tr_step, dt=dt )
#print('CWT-DNN:')
#Y_pred_cwtd = f.pred(modelcwtd, Obs,Obs_labels, tw,tinc, method='cwt',
#                   trace_step=pred_tr_step, dt=dt, freq_index=freq_index)
#print('CWT-CNN:')
#Y_pred_cwtc = f.pred(modelcwtc, Obs,Obs_labels, tw,tinc,method='cwt',
#                   trace_step=pred_tr_step, dt=dt, freq_index=freq_index)

#%
'''plot the result on gather'''
#gather_plot(Obs, tw, savepath='fig/gather_labels_pred.png', 
#            plot_label=True, labels=Y_pred)
fig,ax = f.labels_plot(Y_pred)
ax.set_yticks(np.arange(0,nt/tw+1,400/tw))
ax.set_yticklabels(np.arange(0,nt+1,400)*dt)
#fig.savefig('fig/obs_mlp', dpi=200)

#fig,ax = f.labels_plot(Y_pred_dnn)
fig,ax = f.labels_plot(Y_pred_cnn)
ax.set_yticks(np.arange(0,nt/tw+1,400/tw))
ax.set_yticklabels(np.arange(0,nt+1,400)*dt)
#fig.savefig('fig/obs_cnn', dpi=200)


#%% Results on one trace
'''plot the predicted results on one trace'''

tr_num = 260

obs_samps = f.trace_samp(Obs[:,tr_num], tw, tinc)
f.trace_samp_plot(obs_samps, tw, tinc, nrows, ncols, ann=True, 
                  Y=Obs_labels[:,tr_num], 
                  Y_predit=Y_pred_cnn[:,int(tr_num/pred_tr_step)])

#%% First-break pick -- gather   ################
fb_tr_step = 5
w = 21  # odd, window for attribute calculation
c = 2  # clusters
n = 1  # the neighbouring window for clustering
a = 0  # correction for picked arrival time

'''1-wave classification, 2-clustering inside the first effective wave'''
fb1 = f.fb_pick_gather(Obs, tw, tinc, fb_tr_step, Y_pred_cnn, 
                       pred_tr_step, method='kmeans', w=w, c=c, n=n, a=a)  

#fb2 = f.fb_pick_gather(Obs, tw, tinc, fb_tr_step, Y_pred_cnn, 
#                       pred_tr_step, method='fuzzy', w=w, c=c, n=n) 

'''1-clustering with one whole trace data'''
fb3 = f.fb_pick_gather_wholetrace(Obs, fb_tr_step, method='kmeans', w=w, c=c) 
#fb4 = f.fb_pick_gather_wholetrace(Obs, fb_tr_step, method='fuzzy', w=w, c=c) 


fig, ax = f.seisplot_wig(Obs, inc=fb_tr_step, scale=20, lw=0.2,
                         figsize=(8,6))
ax.scatter(fb1[:,0],fb1[:,1], s=30, 
            facecolors='none',edgecolors='b', lw=1, label='CNN-kmeans')
#ax.scatter(fb2[:,0],fb2[:,1], s=50, marker='D',
#            facecolors='none',edgecolors='b', lw=1, label='fuzzy')
ax.scatter(fb3[:,0],fb3[:,1], s=30, marker='^',
            facecolors='none',edgecolors='r', lw=1, label='Kmeans')
#ax.scatter(fb4[:,0],fb4[:,1], s=60, marker='s',
#            facecolors='none',edgecolors='r', lw=1, label='fuzzy-wholetrace')
ax.legend(loc='lower right' , fontsize=12)
ax.set_yticks(np.arange(0,nt+1,400))
ax.set_yticklabels(np.arange(0,nt+1,400)*dt)
ax.set_xticks(np.arange(0,967,100))
ax.set_xticklabels((np.arange(0,967,100)/fb_tr_step).astype(int))
#fig.savefig('fig/obs_fb', dpi=200)


#%% First-break pick -- one trace

i=260  #13
isamps = f.trace_samp(Obs[:,i], tw, tinc)

'''1-wave classification, 2-clustering inside the first effective wave'''
ifb1, fw1, t1, y1 = f.fb_pick_trace(Obs, i, tw, tinc, Y_pred_cnn, pred_tr_step,
                                    method='kmeans', w=w, c=c, n=n, a=a)
#ifb2, fw2, t2, y2 = f.fb_pick_trace(Obs, i, tw, tinc, Y_pred_cnn, pred_tr_step,
#                                    method='fuzzy', w=w, c=c, n=n)

'''1-clustering using one whole trace data with fuzzy clustering'''
ifb3 = f.fb_pick_wholetrace(Obs[:,i], method='kmeans', w=w, c=c)
#ifb4 = f.fb_pick_wholetrace(Obs[:,i], method='fuzzy', w=31, c=c)
y2 = Obs[:,i][ifb3]
t2 = ifb3%tw

f.trace_samp_fb_plot(isamps, tw,tinc, nrows, ncols, fb_plot=True, fw=fw1, 
                     t1=t1, y1=y1, t2=t2, y2=y2, ann=True, Y=Obs_labels[:,i], 
                     Y_predit=Y_pred_cnn[:,int(i/pred_tr_step)],
                     annx=1, anny=1, ymax=1.5)

fig, ax = plt.subplots(figsize=(13,3))
ax.plot(np.arange(nt), Obs[:,i], 'k',lw=0.8)
ax.scatter(ifb1,0, s=80, facecolors='none',edgecolors='b',lw=2)
#ax.scatter(ifb2,0, s=50, marker='D', facecolors='none',edgecolors='b',lw=2)
ax.scatter(ifb3,0, s=80, marker='o', facecolors='none',edgecolors='r',lw=2)
#ax.scatter(ifb4,0, s=50, marker='s', facecolors='none',edgecolors='y',lw=2)
ax.set_ylim(-0.4,0.4)
ax.set_ylabel('Amplitude', fontsize=13)
ax.set_xlabel('Time samples', fontsize=13)
fig.tight_layout()
ax.legend(loc='upper right' , fontsize=12)
#fig.savefig('fig/s1n_trace_fb', dpi=200)








