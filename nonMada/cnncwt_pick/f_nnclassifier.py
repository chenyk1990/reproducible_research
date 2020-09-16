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
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

import numpy as np
import pandas as pd

from keras.models import Sequential, Model
from keras.layers import Dense, Activation,normalization, Dropout
from keras.layers import Conv1D, Conv2D, Flatten
from keras.layers import AveragePooling2D, MaxPooling1D, MaxPooling2D
from keras.optimizers import SGD, Adam

from keras import callbacks




def dnn(X_train, X_test, Y_train, Y_test, activation, neurons, drop, 
        optimizer, batchsize, epoches, hflayer=2,verbose=1):
    '''
    for DNNs
    the deep feedforward neural network model dealing with 1d vector input
    with only fully connected layers
    '''
    model = Sequential()
    model.add(Dense(units=neurons, input_dim=X_train.shape[1]))
    model.add(Activation(activation)) 
    model.add(Dropout(drop))
    for i in np.arange(hflayer-1):
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))
    model.add(Dense(Y_train.shape[1]))
    model.add(Activation("softmax"))
    
    #binary_accuracy   categorical_crossentropy
    model.compile(loss = 'binary_crossentropy', 
                  optimizer = optimizer, 
                  metrics= ['accuracy'])
    reduce = callbacks.ReduceLROnPlateau(
             monitor='val_loss',factor=0.25, patience=10, mode='min', 
             min_delta=0.005,cooldown=0, min_lr=0.000001)
    earlystop = callbacks.EarlyStopping(monitor='val_loss', 
                                        min_delta=0.001, 
                                        patience=20, 
                                        verbose=0, mode='min')
    print('Training -----------')
    hist = model.fit(X_train, Y_train, 
                     epochs = epoches, 
                     batch_size = batchsize, 
                     verbose=verbose,
                     validation_data = (X_test, Y_test),
                     callbacks = [reduce, earlystop])
    return model, hist


def cnn(X_train, X_test, Y_train, Y_test, activation, neurons, drop, hflayer,
        optimizer, batchsize, epoches, filters=50, kernel_size=5, strides=2,
        pool=True, poolstride=1,verbose=1):
    '''
    for 1D CNNs
    the convolutional neural networks dealing with 1d vector input
    '''
    model = Sequential()
    model.add(Conv2D(input_shape=X_train[0].shape,
                     filters=filters, kernel_size=kernel_size, 
                     strides=strides, padding='same', 
                     data_format='channels_first'))
    if pool==True:
        model.add(AveragePooling2D(pool_size=(1,3), strides=poolstride, 
                                   padding='same',
                               data_format='channels_first'))
        model.add(Dropout(drop))
    
    model.add(Flatten())
    
    for i in np.arange(hflayer):
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop)) 
    
    model.add(Dense(Y_train.shape[1]))
    model.add(Activation("softmax"))
    
    model.compile(loss = 'binary_crossentropy', 
                  optimizer = optimizer, 
                  metrics= ['accuracy'])
    reduce = callbacks.ReduceLROnPlateau(
             monitor='val_loss',factor=0.25, patience=10, mode='min', 
             min_delta=0.005,cooldown=0, min_lr=0.000001)
    earlystop = callbacks.EarlyStopping(monitor='val_loss', 
                                        min_delta=0.001, 
                                        patience=30, 
                                        verbose=0, mode='min')
    print('Training -----------')
    hist = model.fit(X_train, Y_train, 
                     epochs = epoches, 
                     batch_size = batchsize, 
                     verbose=verbose,
                     validation_data = (X_test, Y_test),
                     callbacks = [reduce, earlystop])
    return model, hist


def cwtcnn(X_train, X_test, Y_train, Y_test, 
          activation, neurons, drop, optimizer, batchsize, epoches, hflayer=2,
          reduce_threshold=0.01, lr_patience=10, lr_reduce=0.25, 
          stop_threshold=0.01, stop_patience=25,
          hclayer=1, filters=25, kernel_size=3, strides=2,
          pool=True, poolstride=2,
          verbose=1):
    '''
    the convolutional neural networks dealing with 2d cwt map input
    for CWT-CNNs
    '''
    model = Sequential()
    model.add(Conv2D(input_shape=X_train[0].shape,
                     filters=filters, kernel_size=kernel_size, 
                     strides=strides, padding='same', 
                     data_format='channels_first'))
    for i in np.arange(hclayer-1):
        if pool==True:
            model.add(MaxPooling2D(pool_size=2, strides=poolstride, padding='same',
                                   data_format='channels_first'))
            model.add(Dropout(drop))

        model.add(Conv2D(filters=filters, kernel_size=kernel_size, 
                         strides=strides, padding='same', 
                         data_format='channels_first'))

    #model.add(Activation(activation))
    if pool==True:
        model.add(MaxPooling2D(pool_size=2, strides=poolstride, padding='same',
                               data_format='channels_first'))
    
    model.add(Dropout(drop))
    model.add(Flatten())
    
    for i in np.arange(hflayer-1):
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))       
    
    model.add(Dense(Y_train.shape[1]))
    model.add(Activation("softmax"))
    
    #categorical_crossentropy ,  binary_crossentropy
    model.compile(loss = 'binary_crossentropy', 
                  optimizer = optimizer, 
                  metrics= ['accuracy'])
    reduce = callbacks.ReduceLROnPlateau(
            monitor='val_loss', factor=lr_reduce, patience=lr_patience, 
            mode='min', 
            min_delta=reduce_threshold, cooldown=0, min_lr=0.000001)
    earlystop = callbacks.EarlyStopping(monitor='val_loss', 
                                        min_delta=stop_threshold, 
                                        patience=stop_patience, 
                                        verbose=0, mode='min')
    print('Training -----------')
    hist = model.fit(X_train, Y_train, 
                     epochs = epoches, 
                     batch_size = batchsize, 
                     verbose=verbose,
                     validation_data = (X_test, Y_test),
                     callbacks = [reduce, earlystop])
    return model, hist   


def cwtdnn(X_train, X_test, Y_train, Y_test, 
          activation, neurons, drop, optimizer, batchsize, epoches, hflayer=2,
          reduce_threshold=0.01, lr_patience=10, lr_reduce=0.25, 
          stop_threshold=0.01, stop_patience=25,
          verbose=1):
    '''
    the feedforward neural networks dealing with 2d cwt map input.
    for CWT-DNNs
    '''
    model = Sequential()
    model.add(Flatten(input_shape=X_train[0].shape))        
    
    for i in np.arange(hflayer):
        model.add(Dense(units=neurons))
        model.add(Activation(activation))
        model.add(Dropout(drop))
    
    model.add(Dense(Y_train.shape[1]))
    model.add(Activation("softmax"))
    
    #categorical_crossentropy ,  binary_crossentropy
    model.compile(loss = 'binary_crossentropy', 
                  optimizer = optimizer, 
                  metrics= ['accuracy'])
    reduce = callbacks.ReduceLROnPlateau(
            monitor='val_loss', factor=lr_reduce, patience=lr_patience, 
            mode='min', 
            min_delta=reduce_threshold, cooldown=0, min_lr=0.000001)
    earlystop = callbacks.EarlyStopping(monitor='val_loss', 
                                        min_delta=stop_threshold, 
                                        patience=stop_patience, 
                                        verbose=0, mode='min')
    print('Training -----------')
    hist = model.fit(X_train, Y_train, 
                     epochs = epoches, 
                     batch_size = batchsize, 
                     verbose=verbose,
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
             linewidth=1, c='green', label='Training')
    plt.plot(np.arange(len(cost.val_loss)), cost.val_loss,
             linewidth=1, c='orange', label='Validation')
    plt.xlabel("Epoch",fontsize=14)
    plt.ylabel("Loss",fontsize=14)
    plt.legend(loc='upper right', fontsize=12)
    
    # plotting the accuracy
    plt.subplot(2, 2, 2)
    plt.plot(np.arange(1,len(cost.accuracy)+1), cost.accuracy,
             linewidth=1, c='green', label='Training')
    plt.plot(np.arange(len(cost.val_accuracy)), cost.val_accuracy,
             linewidth=1, c='orange', label='Validating')
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
    plt.show()
