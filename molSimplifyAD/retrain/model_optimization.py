import numpy as np
from functools import partial
from hyperopt import hp, tpe, fmin, Trials, STATUS_OK
from keras.callbacks import EarlyStopping
import tensorflow as tf

from .nets import build_ANN, auc_callback, cal_auc, compile_model


def train_model_hyperopt(hyperspace, X, y, lname,
                         regression=True, epochs=1000,
                         X_val=False, y_val=False, model=False):
    np.random.seed(1234)
    if tf.__version__ >= tf.__version__ >= '2.0.0':
        tf.compat.v1.disable_eager_execution()  ## disable eager in tf2.0 for faster training
    print(("hyperspace: ", hyperspace))
    if not model == False:
        model = build_ANN(hyperspace, X.shape[-1], lname, regression=regression)
    else:
        model = compile_model(model=model, hyperspace=hyperspace,
                              lname=lname, regression=regression)
    if (not X_val == False) and (not y_val == False):
        X_train, X_val = np.split(X, [int(0.8 * X.shape[0])])
        y_train, y_val = np.split(y, [int(0.8 * X.shape[0])])
        y_train, y_val = list(np.transpose(y_train)), list(np.transpose(y_val))
    elif (not X_val == False) or (not y_val == False):
        raise ValueError("Both X_val and y_val need to be specified!")
    else:
        X_train, y_train = X, y
    # if not len(lname) > 1:
    #     y_train, y_val = y_train[0], y_val[0]
    val_data = [X_val, y_val]
    earlystop = EarlyStopping(monitor='val_loss',
                              min_delta=0.0,
                              patience=hyperspace['patience'],
                              verbose=2)
    cb = [earlystop]
    if not regression:
        if len(lname) > 1:
            for ii, ln in enumerate(lname):
                cb += [auc_callback(training_data=(X_train, y_train[ii]), validation_data=(X_val, y_val[ii]), ind=ii)]
        else:
            cb += [auc_callback(training_data=(X_train, y_train), validation_data=(X_val, y_val))]
    history = model.fit(X_train, y_train,
                        epochs=epochs,
                        verbose=2,
                        batch_size=hyperspace['batch_size'],
                        validation_data=val_data,
                        callbacks=cb)
    if regression:
        obj = 0
        if len(lname) > 1:
            for ii, ln in enumerate(lname):
                key = "val_output-%d-%s_scaled_mae" % (ii, ln)
                obj += history.history[key][-1]
            obj /= len(lname)
        else:
            obj = history.history['val_scaled_mae'][-1]
    else:
        if len(lname) > 1:
            val_auc = 0
            for ii, ln in enumerate(lname):
                val_auc += cal_auc(model, X_val, y_val[ii], ind=ii)
        else:
            val_auc = cal_auc(model, X_val, y_val)
        obj = -val_auc / len(lname)
    return {'loss': obj,
            'status': STATUS_OK,
            'epochs': len(history.history[list(history.history.keys())[0]])}


def optimize(X, y, lname,
             regression=True, hyperopt_step=100,
             arch=False, epochs=1000,
             X_val=False, y_val=False,
             model=False):
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
    if not model == False:
        space = {'lr': hp.uniform('lr', 1e-5, 1e-3),
                 'batch_size': hp.choice('batch_size', bzs),
                 'beta_1': hp.uniform('beta_1', 0.75, 0.99),
                 'decay': hp.loguniform('decay', np.log(1e-5), np.log(1e-1)),
                 'amsgrad': True,
                 'patience': 10,
                 }
    else:
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
                             lname=lname,
                             regression=regression,
                             epochs=epochs,
                             X_val=X_val,
                             y_val=y_val,
                             model=model)
    trials = Trials()
    best_params = fmin(objective_func,
                       space,
                       algo=tpe.suggest,
                       trials=trials,
                       max_evals=hyperopt_step,
                       rstate=np.random.RandomState(0)
                       )
    if not model == False:
        best_params.update({'batch_size': bzs[best_params['batch_size']],
                            'amsgrad': True,
                            'patience': 10,
                            })
    else:
        best_params.update({'hidden_size': architectures[best_params['hidden_size']],
                            'batch_size': bzs[best_params['batch_size']],
                            'res': ress[best_params['res']],
                            'bypass': bypasses[best_params['bypass']],
                            'amsgrad': True,
                            'patience': 10,
                            })
    # One extra model training on train/validation set to get the number of epoch for the final model training.
    returned = train_model_hyperopt(best_params, X, y, lname,
                                    regression=regression, epochs=epochs,
                                    X_val=X_val, y_val=y_val,
                                    model=model)
    best_params.update({'epochs': returned['epochs']})
    return best_params
