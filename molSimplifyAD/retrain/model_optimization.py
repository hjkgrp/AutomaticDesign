import numpy as np
from functools import partial
from hyperopt import hp, tpe, fmin, Trials, STATUS_OK
from keras.callbacks import EarlyStopping
import tensorflow as tf

from .nets import build_ANN, auc_callback, cal_auc


def train_model_hyperopt(hyperspace, X, y,
                         regression=True, epochs=1000):
    np.random.seed(1234)
    if tf.__version__ >= tf.__version__ >= '2.0.0':
        tf.compat.v1.disable_eager_execution()  ## disable eager in tf2.0 for faster training
    print("hyperspace: ", hyperspace)
    model = build_ANN(hyperspace, X.shape[-1], regression=regression)
    X_train, X_val = np.split(X, [int(0.8 * X.shape[0])])
    y_train, y_val = np.split(y, [int(0.8 * X.shape[0])])
    val_data = [X_val, y_val]
    earlystop = EarlyStopping(monitor='val_loss',
                              min_delta=0.0,
                              patience=10,
                              verbose=2)
    cb = [earlystop]
    if not regression:
        cb += [auc_callback(training_data=(X_train, y_train), validation_data=(X_val, y_val))]
    history = model.fit(X_train, y_train,
                        epochs=epochs,
                        verbose=2,
                        batch_size=hyperspace['batch_size'],
                        validation_data=val_data,
                        callbacks=cb)
    results = model.evaluate(X_val, y_val)
    if regression:
        obj = results[1]
    else:
        val_auc = cal_auc(model, X_val, y_val)
        obj = -val_auc
    return {'loss': obj,
            'status': STATUS_OK,
            'epochs': len(history.history[list(history.history.keys())[0]])}


def optimize(X, y,
             regression=True, hyperopt_step=100,
             arch=False, epochs=1000):
    np.random.seed(1234)
    if arch == False:
        architectures = [(64,), (128,), (256,), (512,),
                         (64, 64),
                         (128, 128),
                         (256, 256),
                         (512, 512),
                         (64, 64, 64),
                         (128, 128, 128),
                         (256, 256, 256),
                         (512, 512, 512)]
    else:
        architectures = [arch]
    bzs = [16, 32, 64, 128, 256, 512]
    ress = [True, False]
    bypasses = [True, False]
    space = {'lr': hp.uniform('lr', 1e-5, 1e-3),
             'drop_rate': hp.uniform('drop_rate', 0, 0.5),
             'reg': hp.loguniform('reg', np.log(1e-5), np.log(5e-1)),
             'batch_size': hp.choice('batch_size', bzs),
             'hidden_size': hp.choice('hidden_size', architectures),
             'beta_1': hp.uniform('beta_1', 0.75, 0.99),
             'decay': hp.loguniform('decay', np.log(1e-5), np.log(1e-1)),
             'res': hp.choice('res', ress),
             'bypass': hp.choice('bypass', bypasses),
             'amsgrad': True,
             'patience': 10,
             }
    objective_func = partial(train_model_hyperopt,
                             X=X,
                             y=y,
                             regression=regression,
                             epochs=epochs)
    trials = Trials()
    best_params = fmin(objective_func,
                       space,
                       algo=tpe.suggest,
                       trials=trials,
                       max_evals=hyperopt_step,
                       rstate=np.random.RandomState(0)
                       )
    best_params.update({'hidden_size': architectures[best_params['hidden_size']],
                        'batch_size': bzs[best_params['batch_size']],
                        'res': ress[best_params['res']],
                        'bypass': bypasses[best_params['bypass']],
                        'amsgrad': True,
                        'patience': 10,
                        })
    # One extra model training on train/validation set to get the number of epoch for the final model training.
    returned = train_model_hyperopt(best_params, X, y,
                                    regression=regression, epochs=epochs)
    best_params.update({'epochs': returned['epochs']})
    return best_params
