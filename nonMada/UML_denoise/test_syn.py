
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written in Dec, 2019

@author: Yangkang Chen and Mi Zhang
"""
import numpy as np
import tensorflow as tf
from Autoencoder import Autoencoder

def next_batch_random(train_data, batch_size):  
    index = [ i for i in range(0,len(train_data)) ]  
    np.random.shuffle(index);  
    batch_data = []; 
    for i in range(0,batch_size):  
        batch_data.append(train_data[index[i]]);         
    return batch_data

Dtrain=np.fromfile("Datapre/syn_A.dat",dtype=np.float64)
Dtrain=Dtrain.reshape(6000,1600)

Dtest11=np.fromfile("Datapre/syn_A0.dat",dtype=np.float64)
Dtest11=Dtest11.reshape(60,1600)

n_samples=len(Dtrain)
training_epochs =10
batch_size = 100
display_step = 1

corruption_level = 0.0
sparse_reg=0.1
#
n_inputs = 1600
n_hidden = 3200

ae = Autoencoder(n_layers=[n_inputs, n_hidden],
                          transfer_function = tf.nn.sigmoid,
                          optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate = 0.00001),ae_para = [corruption_level, sparse_reg])
init = tf.compat.v1.global_variables_initializer()
sess = tf.compat.v1.Session()
sess.run(init)
for epoch in range(training_epochs):
    avg_cost = 0
    total_batch = int(n_samples / batch_size)
    print("NO of Epochs",training_epochs,"N_samples",n_samples,"NO total_batches",total_batch,"batch size",batch_size)
    # Loop over all batches
    for i in range(total_batch):
#        batch_xs= Dtrain[i*batch_size:(i+1)*batch_size,:] #just consider x, not ylabel
#        ass=mnist.train.images
        batch_xs=next_batch_random(Dtrain,batch_size)
        # Fit training using batch data
        #temp = ae.partial_fit()
        cost, opt = sess.run(ae.partial_fit(), feed_dict={ae.x: batch_xs, ae.keep_prob : ae.in_keep_prob})

        # Compute average loss
        avg_cost += cost / n_samples * batch_size

    # Display logs per epoch step
    if epoch % display_step == 0:
        print("Epoch:", '%d,' % (epoch + 1),
              "Cost:", "{:.9f}".format(avg_cost))
        
ae_test_cost = sess.run(ae.calc_total_cost(), feed_dict={ae.x: Dtest11, ae.keep_prob : ae.in_keep_prob})
print("Total cost: " + str(ae_test_cost))

# Test Process
drecon11 = sess.run(ae.reconstruct(), feed_dict={ae.x: Dtest11, ae.keep_prob : ae.in_keep_prob})
Drecon11=drecon11.astype(np.float)

Drecon11.tofile("Datapre/syn_A1.dat")