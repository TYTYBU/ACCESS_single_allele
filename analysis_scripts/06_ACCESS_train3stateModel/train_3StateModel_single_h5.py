import os
import argparse
import numpy as np
import pandas as pd
import h5py
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, auc, f1_score
import tensorflow as tf
from access_util import *
import logging
import pickle

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = argparse.ArgumentParser(
        description='This script trains 3-state model to predict unbound/bound/recently bound reads.',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    parser.add_argument('h5ad_dir', 
                        type=str, 
                        default=None,
                        help=('Path to h5ad directory containing training and testing data. \n'
                              'Default: None'))

    parser.add_argument('cell_type', 
                        type=str, 
                        default=None,
                        help=('TF cell type. \n'
                              'Default: None'))

    parser.add_argument('TF_name', 
                        type=str, 
                        default=None,
                        help=('TF name. \n'
                              'Default: None'))

    parser.add_argument('motif_str', 
                        type=str, 
                        default=None,
                        help=('TF motif string. \n'
                              'Default: None'))

    parser.add_argument('model_dict_out_dir',
                        type=str,
                        default=None,
                        help=('Output directory for individual model dictionary in pickle format. \n'
                              'Default: None'))

    parser.add_argument('--early_termination',
                        action="store_true", 
                        default=False,
                        help=('Fitting ends early if loss plateaus. \n'
                              'Default: False'))

    parser.add_argument('--num_epochs',
                        type=int,
                        default=30,
                        help=('Number of epochs to fit model. \n'
                              'Default: 30'))

    parser.add_argument('--full_model_only',
                        action='store_true',
                        default=False,
                        help='Only train full model. \n'
                             'Default: False')
    
    parser.add_argument('--no_full_model',
                        action='store_true',
                        default=False,
                        help='Only train base model and edit model. \n'
                             'Default: False')

    return parser.parse_args()

# load arguments
args = parse_args()
h5ad_dir = args.h5ad_dir
cell_type = args.cell_type
TF_name = args.TF_name
motif_str = args.motif_str
model_dict_out_dir = args.model_dict_out_dir

train_base_model = True
train_edit_model = True
train_full_model = True
if args.full_model_only:
    logging.info('Only train full model.')
    train_base_model = False
    train_edit_model = False

if args.no_full_model:
    logging.info('Only train base model and edit model.')
    train_full_model = False
    

BATCH_SIZE = 4096
EPOCHS = args.num_epochs

h5ad_path = h5ad_dir + f'/{cell_type}_{TF_name}_{motif_str}__3readTypes_oneHot.h5'
logging.info(f'Loading h5ad file from {h5ad_path} ...')
hf = h5py.File(h5ad_path, 'r')
x_train = np.array(hf.get('x_train'))
y_train_class = np.array(hf.get('y_train_class'))
x_test = np.array(hf.get('x_test'))
y_test_class = np.array(hf.get('y_test_class'))
hf.close()

logging.info(f'Training data dimensions: {x_train.shape}')
logging.info(f'Testing data dimensions: {x_test.shape}')

# increase input dimension for convolution models
x_train_conv = np.expand_dims(x_train, axis=3) # add one dimension for conv layer
x_test_conv = np.expand_dims(x_test, axis=3) # add one dimension for conv layer

# calculate class weights
samples_per_class = np.array([np.sum(y_train_class==0), np.sum(y_train_class==1), np.sum(y_train_class==2)])
class_weights = sum(samples_per_class) / (len(samples_per_class) * samples_per_class)
class_weights = {i: weight for i, weight in enumerate(class_weights)}
logging.info(f'Class weights: {class_weights}')

# Categorically encode labels
NUM_CLASSES = np.unique(y_train_class).shape[0]
y_train_oneHot = tf.keras.utils.to_categorical(y_train_class, NUM_CLASSES)
y_test_oneHot = tf.keras.utils.to_categorical(y_test_class, NUM_CLASSES)

if train_base_model:
    # Fit data to base model
    logging.info(f'Training base model ...')
    x_train2= x_train_conv[:,0:4,:,:]
    x_test2 = x_test_conv[:,0:4,:,:]
    INPUT_DIM = x_train2.shape[1:]

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=12, restore_best_weights=True)
    lr_scheduler = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=3, min_lr=1e-6)
    if args.early_termination:
        logging.info(f'Early termination is on. Max epochs to run: {args.num_epochs}')
        callback_list = [early_stopping, lr_scheduler]
    else:
        logging.info(f'Early termination is off. Max epochs to run: {args.num_epochs}')
        callback_list = [lr_scheduler]

    model = create_model(input_dim=INPUT_DIM, filter_size=15, output_node_ct=3)
    model.fit(x_train2, y_train_oneHot, batch_size=BATCH_SIZE, epochs=EPOCHS, 
            validation_data=(x_test2, y_test_oneHot), 
            callbacks=callback_list, 
            class_weight=class_weights, verbose=2)
    y_test_pred = model.predict(x_test2, verbose=0)
    y_pred_class = np.argmax(y_test_pred, axis=1)
    f1 = f1_score(y_test_class, y_pred_class, average='weighted')
    logging.info(f'F1 score (weighted): {f1:.3f}')
else:
    model = None


if train_edit_model:
    # Fit data to edit model
    logging.info(f'Training edit model ...')
    x_train2= x_train_conv[:,4:7,:,:]
    x_test2 = x_test_conv[:,4:7,:,:]
    INPUT_DIM = x_train2.shape[1:]

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=12, restore_best_weights=True)
    lr_scheduler = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=3, min_lr=1e-6)
    if args.early_termination:
        logging.info(f'Early termination is on. Max epochs to run: {args.num_epochs}')
        callback_list = [early_stopping, lr_scheduler]
    else:
        logging.info(f'Early termination is off. Max epochs to run: {args.num_epochs}')
        callback_list = [lr_scheduler]

    model2 = create_model(input_dim=INPUT_DIM, filter_size=15, output_node_ct=3)
    model2.fit(x_train2, y_train_oneHot, batch_size=BATCH_SIZE, epochs=EPOCHS, 
            validation_data=(x_test2, y_test_oneHot), 
            callbacks=callback_list, 
            class_weight=class_weights, verbose=2)
    y_test_pred = model2.predict(x_test2, verbose=0)
    y_pred_class = np.argmax(y_test_pred, axis=1)
    f1 = f1_score(y_test_class, y_pred_class, average='weighted')
    logging.info(f'F1 score (weighted): {f1:.3f}')
else:
    model2 = None


if train_full_model:
    # Fit data to full model
    logging.info(f'Training full model ...')
    x_train2= x_train_conv
    x_test2 = x_test_conv
    INPUT_DIM = x_train2.shape[1:]

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=12, restore_best_weights=True)
    lr_scheduler = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=3, min_lr=1e-6)
    if args.early_termination:
        logging.info(f'Early termination is on. Max epochs to run: {args.num_epochs}')
        callback_list = [early_stopping, lr_scheduler]
    else:
        logging.info(f'Early termination is off. Max epochs to run: {args.num_epochs}')
        callback_list = [lr_scheduler]

    model3 = create_model(input_dim=INPUT_DIM, filter_size=15, output_node_ct=3)
    model3.fit(x_train2, y_train_oneHot, batch_size=BATCH_SIZE, epochs=EPOCHS, 
            validation_data=(x_test2, y_test_oneHot), 
            callbacks=callback_list, 
            class_weight=class_weights, verbose=2)
    y_test_pred = model3.predict(x_test2, verbose=0)
    y_pred_class = np.argmax(y_test_pred, axis=1)
    f1 = f1_score(y_test_class, y_pred_class, average='weighted')
    logging.info(f'F1 score (weighted): {f1:.3f}')
else:
    model3 = None

# outputing model dictionary as pickle file
out_path = args.model_dict_out_dir + f'/{cell_type}_{TF_name}_{motif_str}__3stateModels.pkl'
if os.path.exists(out_path):
    logging.info(f'Model file already exists. Overwritting new models ...')
    with open(out_path, 'rb') as file:
        models_dict = pickle.load(file)
    if model is not None:
        models_dict['base_channel'] = model
    if model2 is not None:
        models_dict['edit_channel'] = model2
    if model3 is not None:
        models_dict['all_channel'] = model3

    with open(out_path, 'wb') as handle:
        pickle.dump(models_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    logging.info(f'Model file does not exist. It will be created as {out_path}')
    models_dict = {'base_channel': model, 'edit_channel': model2, 'all_channel': model3}
    with open(out_path, 'wb') as handle:
        pickle.dump(models_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)







