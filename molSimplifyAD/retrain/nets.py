import numpy as np
import tensorflow as tf
from keras import regularizers
from keras.layers import Dense, Dropout, Input, BatchNormalization, Activation, Add, Concatenate
from keras.models import Model
from keras.optimizers import Adam
from sklearn.metrics import roc_auc_score, r2_score
import keras.backend as K
from keras.callbacks import Callback


def scaled_mae(y_true, y_pred):
    return K.mean(K.abs(y_pred - y_true)) / (K.max(y_true) - K.min(y_true))


def mape(y_true, y_pred):
    return 100 * K.mean(K.abs((y_pred - y_true) / (y_true + K.epsilon())))


def r2_val(y_true, y_pred):
    if tf.__version__ >= tf.__version__ >= '2.0.0':
        return tf.py_function(r2_score, (y_true, y_pred), tf.double)
    else:
        return tf.py_func(r2_score, (y_true, y_pred), tf.double)


# def auc(y_true, y_pred):
#     return tf.py_func(roc_auc_score, (y_true, y_pred), tf.double)


def precision(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def recall(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall


def f1(y_true, y_pred):
    p = precision(y_true, y_pred)
    r = recall(y_true, y_pred)
    return 2 * ((p * r) / (p + r + K.epsilon()))


class auc_callback(Callback):
    def __init__(self, training_data, validation_data):
        self.x = training_data[0]
        self.y = training_data[1]
        self.x_val = validation_data[0]
        self.y_val = validation_data[1]

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        y_pred = self.model.predict(self.x)
        roc = roc_auc_score(self.y, y_pred)
        y_pred_val = self.model.predict(self.x_val)
        roc_val = roc_auc_score(self.y_val, y_pred_val)
        print(('roc-auc: %s - roc-auc_val: %s' % (str(round(roc, 4)), str(round(roc_val, 4)))))
        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return


def cal_auc(model, x, y):
    return (roc_auc_score(y, model.predict(x)))


def build_ANN(hyperspace, input_len, regression=True):
    if tf.__version__ >= tf.__version__ >= '2.0.0':
        print("====Tensorflow version >= 2.0.0====")
        model = ANN_tf2(hyperspace, input_len, regression=regression)
    else:
        print("====Tensorflow version < 2.0.0====")
        model = ANN(hyperspace, input_len, regression=regression)
    return model


def ANN(hyperspace, input_len, regression=True):
    np.random.seed(1234)
    inputlayer = Input(shape=(input_len,), name='input')
    layers = [inputlayer]
    for ii in range(len(hyperspace['hidden_size'])):
        layers.append(Dense(hyperspace['hidden_size'][ii],
                            kernel_initializer='he_uniform',
                            kernel_regularizer=regularizers.l2(hyperspace['reg']),
                            name='dense-' + str(ii))(layers[-1]))
        layers.append(BatchNormalization(name='bn-' + str(ii))(layers[-1]))
        layers.append(Activation('relu', name='activation-' + str(ii))(layers[-1]))
        layers.append(Dropout(rate=hyperspace['drop_rate'],
                              name='dropout-' + str(ii))(layers[-1]))
        if hyperspace['res'] and ii:
            base_lookback = -5
            if hyperspace['bypass']:
                base_lookback -= 1
            layers.append(Add(name='sum-' + str(ii))([layers[base_lookback], layers[-1]]))
        if hyperspace['bypass']:
            layers.append(Concatenate(name='concatenate-' + str(ii))([inputlayer, layers[-1]]))
    layers.append(Dense(1, name='dense-last')(layers[-1]))
    outlayer = BatchNormalization(name='bn-last')(layers[-1])
    if not regression:
        outlayer = Activation('sigmoid')(outlayer)
        loss_type = 'binary_crossentropy'
        metrics = ['accuracy', precision, recall, f1]
    else:
        loss_type = 'mse'
        metrics = ['mae', mape, scaled_mae, r2_val]
    model = Model(inputs=[inputlayer], outputs=[outlayer])
    model.compile(loss=loss_type,
                  optimizer=Adam(lr=hyperspace['lr'],
                                 beta_1=hyperspace['beta_1'],
                                 beta_2=0.999,
                                 decay=hyperspace['decay'],
                                 amsgrad=hyperspace['amsgrad']
                                 ),
                  metrics=metrics)
    return model


def ANN_tf2(hyperspace, input_len, regression=True):
    np.random.seed(1234)
    inputlayer = tf.keras.Input(shape=(input_len,), name='input')
    layers = [inputlayer]
    for ii in range(len(hyperspace['hidden_size'])):
        layers.append(tf.keras.layers.Dense(hyperspace['hidden_size'][ii],
                                            kernel_initializer=tf.keras.initializers.he_uniform(),
                                            kernel_regularizer=tf.keras.regularizers.l2(hyperspace['reg']),
                                            name='dense-' + str(ii))(layers[-1]))
        layers.append(tf.keras.layers.BatchNormalization(name='bn-' + str(ii))(layers[-1]))
        layers.append(tf.keras.layers.Activation('relu', name='activation-' + str(ii))(layers[-1]))
        layers.append(tf.keras.layers.Dropout(rate=hyperspace['drop_rate'],
                                              name='dropout-' + str(ii))(layers[-1]))
        if hyperspace['res'] and ii:
            base_lookback = -5
            if hyperspace['bypass']:
                base_lookback -= 1
            layers.append(tf.keras.layers.Add(name='sum-' + str(ii))([layers[base_lookback], layers[-1]]))
        if hyperspace['bypass']:
            layers.append(tf.keras.layers.Concatenate(name='concatenate-' + str(ii))([inputlayer, layers[-1]]))
    layers.append(tf.keras.layers.Dense(1, name='dense-last')(layers[-1]))
    outlayer = tf.keras.layers.BatchNormalization(name='bn-last')(layers[-1])
    if not regression:
        outlayer = tf.keras.layers.Activation('sigmoid', name='activation-last')(outlayer)
        loss_type = tf.keras.losses.BinaryCrossentropy()
        metrics = [tf.keras.metrics.AUC(name="auc"), 'accuracy',
                   tf.keras.metrics.Precision(name="precision"), tf.keras.metrics.Recall(name="recall")]
    else:
        loss_type = 'mse'
        metrics = ['mae', tf.keras.metrics.MeanAbsolutePercentageError(name="mape"),
                   scaled_mae, r2_val]
    model = tf.keras.Model(inputs=[inputlayer], outputs=[outlayer])
    model.compile(loss=loss_type,
                  optimizer=tf.keras.optimizers.Adam(learning_rate=hyperspace['lr'],
                                                     beta_1=hyperspace['beta_1'],
                                                     beta_2=0.999,
                                                     decay=hyperspace['decay'],
                                                     amsgrad=hyperspace['amsgrad']
                                                     ),
                  metrics=metrics)
    return model


### JP's original implementation###
def lr_decay(epoch, initial_lrate=0.1, drop=0.75, epochs_drop=10):
    lrate = initial_lrate * math.pow(drop, math.floor((1 + epoch) / epochs_drop))
    return lrate


def createKerasModel(config_dict, d, p):
    layerlist = []
    fistinput = Input(shape=(d,), name='input')
    layerlist.append(fistinput)
    denseInds = []
    for l in range(1, config_dict['layers']):
        layerlist.append(Dense(config_dict['nodes'], activation=config_dict['activation'], \
                               kernel_regularizer=regularizers.l2(config_dict['l2reg']), \
                               name='dense-' + str(l))(layerlist[-1]))
        layerlist.append(Dropout(config_dict['drop_rate'], name='dropout-' + str(l))(layerlist[-1]))
        if config_dict['semibatch']:
            layerlist.append(BatchNormalization(name='normalization-' + str(l))(layerlist[-1]))

        if config_dict['res'] and not l == 1:
            base_lookback = -3
            if config_dict['semibatch']:
                base_lookback -= 1
            if config_dict['bypass']:
                base_lookback -= 1

            layerlist.append(Add(name='sum-' + str(l))([layerlist[base_lookback], layerlist[-1]]))

        if config_dict['bypass']:
            layerlist.append(Concatenate(name='concatenate-' + str(l))([fistinput, layerlist[-1]]))

            # add output
    layerlist.append(Dense(p, activation=None, name='output')(layerlist[-1]))

    # set model
    model = Model(inputs=layerlist[0], outputs=layerlist[-1])

    # build learning rate scheduler
    this_lr_decay = lambda epoch: lr_decay(epoch, initial_lrate=config_dict["lr"], \
                                           drop=config_dict["decayrate"], \
                                           epochs_drop=config_dict["decayinterval"])

    # build optimizer 
    sgd_optim = keras.optimizers.SGD(lr=this_lr_decay(0), momentum=config_dict['momentum'],
                                     decay=0.0, nesterov=config_dict['nest'])

    # compile model
    earlystop = EarlyStopping(monitor='val_mean_absolute_error', min_delta=config_dict['min_delta'],
                              patience=config_dict['patience'], mode='auto', verbose=1)

    # compile callbacks
    callback_list = [LearningRateScheduler(this_lr_decay, verbose=0), earlystop]

    # compile model
    model.compile(loss='mse', optimizer=sgd_optim, metrics=['mse', 'mae'])
    # model.compile(loss=quarticloss, optimizer=sgd_optim,metrics=['mse', 'mae'])

    return (model, callback_list)
