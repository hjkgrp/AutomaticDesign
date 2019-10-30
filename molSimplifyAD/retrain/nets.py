import numpy as np
from keras import regularizers
from keras.layers import Dense, Dropout, Input, BatchNormalization, Activation, Add, Concatenate
from keras.models import Model
from keras.optimizers import Adam


def ANN(hyperspace, input_len, regression=True):
    np.random.seed(1234)
    inputlayer = Input(shape=(input_len,), name='input')
    layers = [inputlayer]
    for ii in range(len(hyperspace['hidden_size'])):
        layers.append(Dense(hyperspace['hidden_size'][ii],
                            kernel_initializer='he_uniform',
                            kernel_regularizer=regularizers.l2(hyperspace['reg']),
                            name='dense-'+str(ii))(layers[-1]))
        layers.append(BatchNormalization(name='bn-'+str(ii))(layers[-1]))
        layers.append(Activation('relu', name='activation-'+str(ii))(layers[-1]))
        layers.append(Dropout(rate=hyperspace['drop_rate'],
                              name='dropout-'+str(ii))(layers[-1]))
        if hyperspace['res'] and ii:
            base_lookback = -5
            if hyperspace['bypass']:
                base_lookback -= 1
            layers.append(Add(name='sum-'+str(ii))([layers[base_lookback], layers[-1]]))
        if hyperspace['bypass']:
            layers.append(Concatenate(name='concatenate-'+str(ii))([inputlayer, layers[-1]]))
    layers.append(Dense(1, name='dense-last')(layers[-1]))
    outlayer = BatchNormalization(name='bn-last')(layers[-1])
    if not regression:
        outlayer = Activation('sigmoid')(outlayer)
        loss_type = 'binary_crossentropy'
        metrics = ['accuracy']
    else:
        loss_type = 'mse'
        metrics = ['mae']
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


### JP's original implementation###
def lr_decay(epoch,initial_lrate = 0.1, drop=0.75, epochs_drop=10):
    lrate = initial_lrate * math.pow(drop, math.floor((1+epoch)/epochs_drop))
    return lrate

def createKerasModel(config_dict,d,p):
    layerlist= []
    fistinput = Input(shape=(d, ), name='input')
    layerlist.append(fistinput)
    denseInds = []
    for l in range(1,config_dict['layers']):
        layerlist.append(Dense(config_dict['nodes'],activation=config_dict['activation'],\
                               kernel_regularizer=regularizers.l2(config_dict['l2reg']),\
                                name='dense-'+str(l))(layerlist[-1]))
        layerlist.append(Dropout(config_dict['drop_rate'],name='dropout-'+str(l))(layerlist[-1]))
        if config_dict['semibatch']:
            layerlist.append(BatchNormalization(name='normalization-'+str(l))(layerlist[-1]))



        if config_dict['res'] and not l==1:
            base_lookback = -3
            if config_dict['semibatch']:
                base_lookback -= 1
            if config_dict['bypass']:
                base_lookback -= 1

            layerlist.append(Add(name='sum-'+str(l))([layerlist[base_lookback], layerlist[-1]]))        

        if config_dict['bypass']:
            layerlist.append(Concatenate(name='concatenate-'+str(l))([fistinput, layerlist[-1]]))                   

    # add output 
    layerlist.append(Dense(p,activation=None, name='output')(layerlist[-1]))

    # set model
    model = Model(inputs=layerlist[0], outputs=layerlist[-1])


    # build learning rate scheduler
    this_lr_decay = lambda epoch: lr_decay(epoch,initial_lrate = config_dict["lr"],\
                                           drop=config_dict["decayrate"], \
                                           epochs_drop=config_dict["decayinterval"])

    # build optimizer 
    sgd_optim =  keras.optimizers.SGD(lr=this_lr_decay(0), momentum=config_dict['momentum'],
                                      decay=0.0, nesterov=config_dict['nest'])




    # compile model
    earlystop = EarlyStopping(monitor='val_mean_absolute_error', min_delta=config_dict['min_delta'],
                              patience=config_dict['patience'],mode='auto',verbose=1)

    # compile callbacks
    callback_list = [LearningRateScheduler(this_lr_decay,verbose=0),earlystop]

    # compile model
    model.compile(loss='mse', optimizer=sgd_optim,metrics=['mse', 'mae'])
    #model.compile(loss=quarticloss, optimizer=sgd_optim,metrics=['mse', 'mae'])
    
    return(model, callback_list)