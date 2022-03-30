import numpy as np
import tensorflow as tf
# from tensorflow.keras import layers
from tensorflow.keras.models import Model
from tensorflow.keras import layers, losses
import matplotlib.pyplot as plt
# tf.compat.v1.disable_eager_execution()

latent_dim = 128 

class Autoencoder(Model):
  def __init__(self, latent_dim):
    super(Autoencoder, self).__init__()
    self.latent_dim = latent_dim   
    self.encoder = tf.keras.Sequential([
      layers.Dense(latent_dim, activation='relu'),
    ])
    self.decoder = tf.keras.Sequential([
      layers.Dense(40*40, activation='sigmoid'),
    ])

  def call(self, x):
    encoded = self.encoder(x)
    decoded = self.decoder(encoded)
    return decoded


def display(array1, array2):
    """
    Displays ten random images from each one of the supplied arrays.
    """

    n = 10

    indices = np.random.randint(len(array1), size=n)
    images1 = array1[indices, :]
    images2 = array2[indices, :]

    plt.figure(figsize=(20, 4))
    for i, (image1, image2) in enumerate(zip(images1, images2)):
        ax = plt.subplot(2, n, i + 1)
        plt.imshow(image1.reshape(40, 40))
        plt.gray()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        ax = plt.subplot(2, n, i + 1 + n)
        plt.imshow(image2.reshape(40, 40))
        plt.gray()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    plt.show()
    
Dtrain=np.fromfile("Datapre/syn_A.dat",dtype=np.float64)
Dtrain=Dtrain.reshape(6000,40*40)

Dtest11=np.fromfile("Datapre/syn_A0.dat",dtype=np.float64)
Dtest11=Dtest11.reshape(60,40*40)

autoencoder = Autoencoder(latent_dim)
autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())

autoencoder.fit(Dtrain, Dtrain,
                epochs=20,
                shuffle=True,
                validation_data=(Dtest11, Dtest11))

encoded_imgs = autoencoder.encoder(Dtrain).numpy()
decoded_imgs = autoencoder.decoder(encoded_imgs).numpy()

Drecon=decoded_imgs.astype(np.float)
Drecon.tofile("Datapre/syn_A1.dat")

display(Dtrain, Drecon)
                
                