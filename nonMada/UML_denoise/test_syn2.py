import numpy as np
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.models import Model
import matplotlib.pyplot as plt
# tf.compat.v1.disable_eager_execution()

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
Dtrain=Dtrain.reshape(6000,40,40,1)

Dtest11=np.fromfile("Datapre/syn_A0.dat",dtype=np.float64)
Dtest11=Dtest11.reshape(60,40,40,1)


input = layers.Input(shape=(40, 40, 1))

# Encoder
x = layers.Conv2D(32, (3, 3), activation="relu", padding="same")(input)
x = layers.MaxPooling2D((2, 2), padding="same")(x)
x = layers.Conv2D(32, (3, 3), activation="relu", padding="same")(x)
x = layers.MaxPooling2D((2, 2), padding="same")(x)

# Decoder
x = layers.Conv2DTranspose(32, (3, 3), strides=2, activation="relu", padding="same")(x)
x = layers.Conv2DTranspose(32, (3, 3), strides=2, activation="relu", padding="same")(x)
x = layers.Conv2D(1, (3, 3), activation="sigmoid", padding="same")(x)

# Autoencoder
autoencoder = Model(input, x)
autoencoder.compile(optimizer="adam", loss="binary_crossentropy")
autoencoder.summary()

autoencoder.fit(
    x=Dtrain,
    y=Dtrain,
    epochs=50,
    batch_size=128,
    shuffle=True,
    validation_data=(Dtest11, Dtest11),
)


predictions = autoencoder.predict(Dtest11)
display(Dtest11, predictions)

drecon= autoencoder.predict(Dtrain)
Drecon=drecon.astype(np.float)
Drecon.tofile("Datapre/syn_A1.dat")


display(Dtrain, Drecon)

