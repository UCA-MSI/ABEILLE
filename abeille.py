# -*- coding: utf-8 -*-
#!/usr/bin/python

import time
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
tfk = tf.keras
tfkl = tfk.layers


def abeille_VAE(file, logarithm = True, read_count = True,
                batch_size = 30, epochs = 1500, kernel = "lecun_normal",
                kl_loss_weight = 0.5, edl1 = 2048, edl2 = 1024, edl3 = 512,
                edl4 = 256, latent_size = 128, ddl1 = 256, ddl2 = 512,
                ddl3 = 1024, ddl4 = 2048):
    def sampling(args):
        z_mean, z_log_var = args
        batch = tfk.backend.shape(z_mean)[0]
        dim = tfk.backend.int_shape(z_mean)[1]
        epsilon = tfk.backend.random_normal(shape=(batch, dim))
        out = z_mean + tf.keras.backend.exp(0.5 * z_log_var) * epsilon
        return out
    data = pd.read_csv(file, index_col = 0)
    if logarithm == True:
        data_input = np.log(data.T+1)
    else:
        data_input = data.T
    
    n_cols = data_input.shape[1]    
    start_time = time.time()
    #####ENCODER
    #Input
    inputs = tfk.Input(shape=(n_cols,), name='encoder_input')
    #First hidden layer
    hidden_encoder_1 = tfkl.Dense(edl1, kernel_initializer = kernel, name="hidden_encoder_1")(inputs)
    encoder_norm_1 = tfkl.BatchNormalization(name="encoder_norm_1")(hidden_encoder_1)
    hidden_encoder_1_activation = tfkl.ELU()(encoder_norm_1)
    #Second hidden layer
    hidden_encoder_2 = tfkl.Dense(edl2, kernel_initializer = kernel, name="hidden_encoder_2")(hidden_encoder_1_activation)
    encoder_norm_2 = tfkl.BatchNormalization(name="encoder_norm_2")(hidden_encoder_2)
    hidden_encoder_2_activation = tfkl.ELU()(encoder_norm_2)
    #Third hidden layer
    hidden_encoder_3 = tfkl.Dense(edl3, kernel_initializer = kernel, name="hidden_encoder_3")(hidden_encoder_2_activation)
    encoder_norm_3 = tfkl.BatchNormalization(name="encoder_norm_3")(hidden_encoder_3)
    hidden_encoder_3_activation = tfkl.ELU()(encoder_norm_3)
    #Fourth hidden layer
    hidden_encoder_4 = tfkl.Dense(edl4, kernel_initializer = kernel, name="hidden_encoder_4")(hidden_encoder_3_activation)
    encoder_norm_4 = tfkl.BatchNormalization(name="encoder_norm_4")(hidden_encoder_4)
    hidden_encoder_4_activation = tfkl.ELU()(encoder_norm_4)
    #Mean for the sampling
    z_mean = tfkl.Dense(latent_size, name='z_mean')(hidden_encoder_4_activation)
    #Var for the sampling
    z_log_var = tfkl.Dense(latent_size, name='z_log_var')(hidden_encoder_4_activation)
    #Sample from the values below
    z = tfkl.Lambda(sampling, output_shape=(latent_size,), name='z')([z_mean, z_log_var])
    #Encoder model
    encoder = tfk.Model(inputs, [z_mean, z_log_var, z], name='encoder')
    
    #####DECODER
    #Input
    latent_inputs = tfk.Input(shape=(latent_size,), name='z_sampling')
    #Fifth hidden layer
    hidden_decoder_5 = tfkl.Dense(ddl1, kernel_initializer = kernel, name="hidden_decoder_4")(latent_inputs)
    decoder_norm_5 = tfkl.BatchNormalization(name="decoder_norm_4")(hidden_decoder_5)
    hidden_encoder_5_activation = tfkl.ELU()(decoder_norm_5)
    #Sixth hidden layer
    hidden_decoder_6 = tfkl.Dense(ddl2, kernel_initializer = kernel, name="hidden_decoder_5")(hidden_encoder_5_activation)
    decoder_norm_6 = tfkl.BatchNormalization(name="decoder_norm_5")(hidden_decoder_6)
    hidden_encoder_6_activation = tfkl.ELU()(decoder_norm_6)
    #Seventh hidden layer
    hidden_decoder_7 = tfkl.Dense(ddl3, kernel_initializer = kernel, name="hidden_decoder_6")(hidden_encoder_6_activation)
    decoder_norm_7 = tfkl.BatchNormalization(name="decoder_norm_6")(hidden_decoder_7)
    hidden_encoder_7_activation = tfkl.ELU()(decoder_norm_7)
    #Eighth hidden layer
    hidden_decoder_8 = tfkl.Dense(ddl4, kernel_initializer = kernel, name="hidden_decoder_7")(hidden_encoder_7_activation)
    decoder_norm_8 = tfkl.BatchNormalization(name="decoder_norm_7")(hidden_decoder_8)
    hidden_encoder_8_activation = tfkl.ELU()(decoder_norm_8)
    #Output
    outputs = tfkl.Dense(n_cols)(hidden_encoder_8_activation)
    #Decoder model
    decoder = tfk.Model(latent_inputs, outputs, name="decoder")
    
    #####VAE
    outputs = decoder(encoder(inputs)[2])
    vae = tfk.Model(inputs, outputs, name="vae")
    
    #####Personalized loss
    msle = tfk.losses.MSLE(inputs, outputs)
    msle *= n_cols
    kl_loss = 1 + z_log_var - tfk.backend.square(z_mean) - tfk.backend.exp(z_log_var)
    kl_loss = tfk.backend.sum(kl_loss, axis=-1)
    kl_loss *= -kl_loss_weight
    vae_loss = tfk.backend.mean(msle + kl_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer="adam")
    vae.summary()
    
    #####Train
    vae.fit(data_input, epochs=epochs, batch_size=batch_size)
    end_time = time.time() - start_time
    print("-----",end_time//3600,"hours",(end_time%3600)//60,"minutes",round((end_time%3600)%60,2),"secondes","-----")
    
    pred = vae.predict(data_input)
    pred = pd.DataFrame(pred)
    pred.index = data_input.index
    pred.columns = data_input.columns
    if logarithm == True:
        pred = np.exp(pred) - 1
    if read_count == True:
        pred = pred * (pred > 0)
    pred = pred.T
    return(pred)


if __name__ == "__main__":

    args = parser = argparse.ArgumentParser(
        description='Launch ABEILLE VAE')
    parser.add_argument('file', type=str, help='File to reconstruct')
    parser.add_argument('output', type = str, help = 'Output file name')

    args = parser.parse_args()

    data_recons = abeille_VAE(args.file)
    data_recons.to_csv(args.output)