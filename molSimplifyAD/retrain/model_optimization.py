import numpy as np
from functools import partial
from hyperopt import hp, tpe, fmin, Trials, STATUS_OK
from keras.callbacks import EarlyStopping

from nets import ANN


def train_model_hyperopt(hyperspace, X, y, regression=True):
    np.random.seed(1234)
    print("hyperspace: ", hyperspace)
    model = ANN(hyperspace, X.shape[-1], regression=regression)
    X_train, X_val = np.split(X, [int(0.8 * X.shape[0])])
    y_train, y_val = np.split(y, [int(0.8 * X.shape[0])])
    val_data = [X_val, y_val]
    earlystop = EarlyStopping(monitor='val_loss',
                              min_delta=0.0,
                              patience=10,
                              verbose=2)
    cb = [earlystop]
    history = model.fit(X_train, y_train,
                        epochs=1000,
                        verbose=2,
                        batch_size=hyperspace['batch_size'],
                        validation_data=val_data,
                        callbacks=cb)
    results = model.evaluate(X_val, y_val)
    if regression:
        obj = results[1]
    else:
        obj = -results[1]
    return {'loss': obj,
            'status': STATUS_OK,
            'epochs': len(history.history[history.history.keys()[0]])}


def optimize(X, y, regression=True):
    np.random.seed(1234)
    architectures = [(64,), (128,), (256,), (512,),
                     (64, 64),
                     (128, 128),
                     (256, 256),
                     (512, 512),
                     (64, 64, 64),
                     (128, 128, 128),
                     (256, 256, 256),
                     (512, 512, 512)]
    bzs = [16, 32, 64, 128, 256, 512]
    space = {'lr': hp.uniform('lr', 1e-5, 1e-3),
             'drop_rate': hp.uniform('drop_rate', 0, 0.5),
             'reg': hp.loguniform('reg', np.log(1e-5), np.log(5e-1)),
             'batch_size': hp.choice('batch_size', bzs),
             'hidden_size': hp.choice('hidden_size', architectures),
             'beta_1': hp.uniform('beta_1', 0.75, 0.99),
             'decay': hp.loguniform('decay', np.log(1e-5), np.log(1e-1)),
             'amsgrad': True,
             'patience': 10,
             }
    objective_func = partial(train_model_hyperopt,
                             X=X,
                             y=y,
                             regression=regression)
    trials = Trials()
    best_params = fmin(objective_func,
                       space,
                       algo=tpe.suggest,
                       trials=trials,
                       max_evals=200,
                       rstate=np.random.RandomState(0)
                       )
    best_params.update({'hidden_size': architectures[best_params['hidden_size']],
                        'batch_size': bzs[best_params['batch_size']],
                        'amsgrad': True,
                        'patience': 10,
                        })
    returned = train_model_hyperopt(best_params, X, y, regression=regression)
    best_params.update({'epochs': returned['epochs']})
    return best_params
