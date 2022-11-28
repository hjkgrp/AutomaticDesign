import numpy as np
import tensorflow as tf
from keras import regularizers
from keras.layers import Dense, Dropout, Input, BatchNormalization, Activation, Add, Concatenate, LeakyReLU
from keras.models import Model
from keras.optimizers import Adam
from sklearn.metrics import roc_auc_score, r2_score
import keras.backend as K
from keras.callbacks import Callback


def scaled_mae(y_true, y_pred):
    return K.mean(K.abs(y_pred - y_true)) / (K.max(y_true) - K.min(y_true) + K.epsilon())


def mape(y_true, y_pred):
    return 100 * K.mean(K.abs((y_pred - y_true) / (y_true + K.epsilon())))


def r2_val(y_true, y_pred):
    if tf.__version__ >= tf.__version__ >= '2.0.0':
        return tf.py_function(r2_score, (y_true, y_pred), tf.double)
    else:
        return tf.py_func(r2_score, (y_true, y_pred), tf.double)


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
    ## Only for tensorflow 1.x where the AUC metric is not available.
    def __init__(self, training_data, validation_data, ind=None):
        self.x = training_data[0]
        self.y = np.array(training_data[1])
        self.x_val = validation_data[0]
        self.y_val = np.array(validation_data[1])
        self.ind = ind

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        if not self.ind == None:
            y_pred = self.model.predict(self.x)[self.ind]
            y_pred_val = self.model.predict(self.x_val)[self.ind]
            ii = self.ind
        else:
            y_pred = self.model.predict(self.x)
            y_pred_val = self.model.predict(self.x_val)
            ii = 0
        # print("111: ", self.y.reshape(-1, 1).shape, y_pred.shape)
        roc = roc_auc_score(self.y.reshape(-1, 1), y_pred)
        roc_val = roc_auc_score(self.y_val.reshape(-1, 1), y_pred_val)
        print(('%d: roc-auc: %s - roc-auc_val: %s' % (ii, str(round(roc, 4)), str(round(roc_val, 4)))))
        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return


def cal_auc(model, x, y, ind=None):
    if not ind == None:
        return roc_auc_score(np.array(y).reshape(-1, 1), model.predict(x)[ind])
    else:
        print(np.array(y).shape)
        return roc_auc_score(np.array(y).reshape(-1, 1), model.predict(x))


def build_ANN(hyperspace, input_len, lname, regression=True):
    '''
    Build a fully-connected neural network.

    Parameters
    ---
        hyperspace: dict, a dictionary of hyperparameters
        input_len: int, length of the feeature vector
        lname: str, name of the target property
        regression: boolean, whether it is a regression task

    Returns
    ---
        model: keras.model object
    '''
    if tf.__version__ >= '2.0.0':
        print("====Tensorflow version >= 2.0.0====")
        model = ANN_tf2(hyperspace, input_len, lname, regression=regression)
    else:
        print("====Tensorflow version < 2.0.0====")
        model = ANN(hyperspace, input_len, lname, regression=regression)
    return model


def ANN(hyperspace, input_len, lname, regression=True):
    np.random.seed(1234)
    inputlayer = Input(shape=(input_len,), name='input')
    layers = [inputlayer]
    for ii in range(len(hyperspace['hidden_size'])):
        layers.append(Dense(hyperspace['hidden_size'][ii],
                            kernel_initializer='he_uniform',
                            kernel_regularizer=regularizers.l2(hyperspace['reg']),
                            name='dense-' + str(ii))(layers[-1]))
        layers.append(BatchNormalization(name='bn-' + str(ii))(layers[-1]))
        layers.append(LeakyReLU(alpha=0.2, name='activation-' + str(ii))(layers[-1]))
        layers.append(Dropout(rate=hyperspace['drop_rate'],
                              name='dropout-' + str(ii))(layers[-1]))
        if ('res' in hyperspace.keys()) and (hyperspace['res'] and ii):
            base_lookback = -5
            if hyperspace['bypass']:
                base_lookback -= 1
            layers.append(Add(name='sum-' + str(ii))([layers[base_lookback], layers[-1]]))
        if ('bypass' in hyperspace.keys()) and hyperspace['bypass']:
            layers.append(Concatenate(name='concatenate-' + str(ii))([inputlayer, layers[-1]]))
    last_dense, outlayer, loss_weights, loss_type = [], [], [], []
    for ii, ln in enumerate(lname):
        last_dense.append(Dense(1, name='dense-last-%d' % ii)(layers[-1]))
        outlayer.append(0)
        if not regression:
            last_bn = BatchNormalization(name='bn-last-%d' % ii)(last_dense[ii])
            outlayer[ii] = Activation('sigmoid', name='output-%d-%s' % (ii, ln))(last_bn)
            _loss_type = 'binary_crossentropy'
            metrics = ['accuracy', precision, recall, f1]
        else:
            outlayer[ii] = BatchNormalization(name='output-%d-%s' % (ii, ln))(last_dense[ii])
            _loss_type = 'mse'
            metrics = ['mae', scaled_mae, r2_val]
        loss_weights.append(1.0)
        loss_type.append(_loss_type)
    model = Model(inputs=[inputlayer], outputs=outlayer)
    model.compile(loss=loss_type,
                  optimizer=Adam(lr=hyperspace['lr'],
                                 beta_1=hyperspace['beta_1'],
                                 beta_2=0.999,
                                 decay=hyperspace['decay'],
                                 amsgrad=hyperspace['amsgrad']
                                 ),
                  metrics=metrics,
                  loss_weights=loss_weights,
                  )
    return model


def ANN_tf2(hyperspace, input_len, lname, regression=True):
    np.random.seed(1234)
    inputlayer = tf.keras.Input(shape=(input_len,), name='input')
    layers = [inputlayer]
    for ii in range(len(hyperspace['hidden_size'])):
        layers.append(tf.keras.layers.Dense(hyperspace['hidden_size'][ii],
                                            kernel_initializer=tf.keras.initializers.he_uniform(),
                                            kernel_regularizer=tf.keras.regularizers.l2(hyperspace['reg']),
                                            name='dense-' + str(ii))(layers[-1]))
        layers.append(tf.keras.layers.BatchNormalization(name='bn-' + str(ii))(layers[-1]))
        layers.append(tf.keras.layers.LeakyReLU(alpha=0.2, name='activation-' + str(ii))(layers[-1]))
        layers.append(tf.keras.layers.Dropout(rate=hyperspace['drop_rate'],
                                              name='dropout-' + str(ii))(layers[-1]))
        if hyperspace['res'] and ii:
            base_lookback = -5
            if hyperspace['bypass']:
                base_lookback -= 1
            layers.append(tf.keras.layers.Add(name='sum-' + str(ii))([layers[base_lookback], layers[-1]]))
        if hyperspace['bypass']:
            layers.append(tf.keras.layers.Concatenate(name='concatenate-' + str(ii))([inputlayer, layers[-1]]))
    last_dense, outlayer, loss_weights, loss_type = [], [], [], []
    for ii, ln in enumerate(lname):
        last_dense.append(Dense(1, name='dense-last-%d' % ii)(layers[-1]))
        outlayer.append(0)
        if not regression:
            last_bn = BatchNormalization(name='bn-last-%d' % ii)(last_dense[ii])
            outlayer[ii] = tf.keras.layers.Activation('sigmoid', name='activation-last')(last_bn)
            _loss_type = tf.keras.losses.BinaryCrossentropy()
            metrics = [tf.keras.metrics.AUC(name="auc"), 'accuracy',
                       tf.keras.metrics.Precision(name="precision"), tf.keras.metrics.Recall(name="recall")]
        else:
            outlayer[ii] = BatchNormalization(name='output-%d-%s' % (ii, ln))(last_dense[ii])
            _loss_type = 'mse'
            metrics = ['mae', tf.keras.metrics.MeanAbsolutePercentageError(name="mape"), 
                       scaled_mae, r2_val]
        loss_weights.append(1.0)
        loss_type.append(_loss_type)
    model = tf.keras.Model(inputs=[inputlayer], outputs=outlayer)
    model.compile(loss=loss_type,
                  optimizer=tf.keras.optimizers.Adam(learning_rate=hyperspace['lr'],
                                                     beta_1=hyperspace['beta_1'],
                                                     beta_2=0.999,
                                                     decay=hyperspace['decay'],
                                                     amsgrad=hyperspace['amsgrad']
                                                     ),
                  metrics=metrics,
                  loss_weights=loss_weights,)
    return model


def compile_model(model, hyperspace, lname, regression):
    last_dense, outlayer, loss_weights, loss_type = [], [], [], []
    for ii, ln in enumerate(lname):
        if not regression:
            _loss_type = 'binary_crossentropy'
            metrics = ['accuracy', precision, recall, f1]
        else:
            _loss_type = 'mean_squared_error'
            metrics = ['mean_absolute_error', scaled_mae, r2_val]
        loss_weights.append(1.0)
        loss_type.append(_loss_type)
    model.compile(loss=loss_type,
                  optimizer=Adam(lr=hyperspace['lr'],
                                 beta_1=hyperspace['beta_1'],
                                 beta_2=0.999,
                                 decay=hyperspace['decay'],
                                 amsgrad=hyperspace['amsgrad']
                                 ),
                  metrics=metrics,
                  loss_weights=loss_weights,
                  )
    return model


### JP's original implementation###
def lr_decay(epoch, initial_lrate=0.1, drop=0.75, epochs_drop=10):
    lrate = initial_lrate * math.pow(drop, math.floor((1 + epoch) / epochs_drop))
    return lrate


def createKerasModel(config_dict, d, p):
    ### JP's old function, can probably be deprecated.
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
