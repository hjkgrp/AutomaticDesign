import numpy as np
from keras import regularizers
from keras.layers import Dense, Dropout, Input, BatchNormalization, Activation
from keras.models import Model
from keras.optimizers import Adam


def ANN(hyperspace, input_len, regression=True):
    np.random.seed(1234)
    input = Input(shape=(input_len,), name='input')
    m = Dense(hyperspace['hidden_size'][0],
              kernel_initializer='he_uniform',
              kernel_regularizer=regularizers.l2(hyperspace['reg']))(input)
    m = BatchNormalization()(m)
    m = Activation('relu')(m)
    m = Dropout(rate=hyperspace['drop_rate'])(m)
    for ii in range(1, len(hyperspace['hidden_size'])):
        m = Dense(hyperspace['hidden_size'][ii],
                  kernel_initializer='he_uniform',
                  kernel_regularizer=regularizers.l2(hyperspace['reg']))(m)
        m = BatchNormalization()(m)
        m = Activation('relu')(m)
        m = Dropout(rate=hyperspace['drop_rate'])(m)
    m = Dense(1)(m)
    output = BatchNormalization()(m)
    if not regression:
        output = Activation('sigmoid')(output)
        loss_type = 'binary_crossentropy'
        metrics = ['accuracy']
    else:
        loss_type = 'mse'
        metrics = ['mae']
    model = Model(inputs=[input, ], outputs=[output])
    model.compile(loss=loss_type,
                  optimizer=Adam(lr=hyperspace['lr'],
                                 beta_1=hyperspace['beta_1'],
                                 beta_2=0.999,
                                 decay=hyperspace['decay'],
                                 amsgrad=hyperspace['amsgrad']
                                 ),
                  metrics=metrics)
    return model
