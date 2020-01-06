import numpy as np
import pandas as pd
import tensorflow as tf
from datetime import datetime
from keras.models import load_model
import sklearn.preprocessing
import sklearn.utils
from molSimplifyAD.utils.pymongo_tools import convert2dataframe, connect2db, push_models
from pkg_resources import resource_filename, Requirement
from molSimplify.python_nn.tf_ANN import get_key, load_ANN_variables, load_keras_ann, initialize_model_weights
from .nets import build_ANN, cal_auc
from .model_optimization import optimize
from .calculate_pairing_prop import group_conditions, pairing

name_converter_dict = {"oxstate": "ox", "spinmult": "spin", "charge_lig": "ligcharge"}
RACs180 = ['D_lc-I-0-ax', 'D_lc-I-0-eq', 'D_lc-I-1-ax', 'D_lc-I-1-eq', 'D_lc-I-2-ax', 'D_lc-I-2-eq', 'D_lc-I-3-ax',
           'D_lc-I-3-eq', 'D_lc-S-0-ax', 'D_lc-S-0-eq', 'D_lc-S-1-ax', 'D_lc-S-1-eq', 'D_lc-S-2-ax', 'D_lc-S-2-eq',
           'D_lc-S-3-ax', 'D_lc-S-3-eq', 'D_lc-T-0-ax', 'D_lc-T-0-eq', 'D_lc-T-1-ax', 'D_lc-T-1-eq', 'D_lc-T-2-ax',
           'D_lc-T-2-eq', 'D_lc-T-3-ax', 'D_lc-T-3-eq', 'D_lc-Z-0-ax', 'D_lc-Z-0-eq', 'D_lc-Z-1-ax', 'D_lc-Z-1-eq',
           'D_lc-Z-2-ax', 'D_lc-Z-2-eq', 'D_lc-Z-3-ax', 'D_lc-Z-3-eq', 'D_lc-chi-0-ax', 'D_lc-chi-0-eq',
           'D_lc-chi-1-ax', 'D_lc-chi-1-eq', 'D_lc-chi-2-ax', 'D_lc-chi-2-eq', 'D_lc-chi-3-ax', 'D_lc-chi-3-eq',
           'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all', 'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all',
           'D_mc-S-2-all', 'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all', 'D_mc-T-3-all',
           'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all', 'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all',
           'D_mc-chi-2-all', 'D_mc-chi-3-all', 'f-I-0-all', 'f-I-0-ax', 'f-I-0-eq', 'f-I-1-all', 'f-I-1-ax', 'f-I-1-eq',
           'f-I-2-all', 'f-I-2-ax', 'f-I-2-eq', 'f-I-3-all', 'f-I-3-ax', 'f-I-3-eq', 'f-S-0-all', 'f-S-0-ax',
           'f-S-0-eq', 'f-S-1-all', 'f-S-1-ax', 'f-S-1-eq', 'f-S-2-all', 'f-S-2-ax', 'f-S-2-eq', 'f-S-3-all',
           'f-S-3-ax', 'f-S-3-eq', 'f-T-0-all', 'f-T-0-ax', 'f-T-0-eq', 'f-T-1-all', 'f-T-1-ax', 'f-T-1-eq',
           'f-T-2-all', 'f-T-2-ax', 'f-T-2-eq', 'f-T-3-all', 'f-T-3-ax', 'f-T-3-eq', 'f-Z-0-all', 'f-Z-0-ax',
           'f-Z-0-eq', 'f-Z-1-all', 'f-Z-1-ax', 'f-Z-1-eq', 'f-Z-2-all', 'f-Z-2-ax', 'f-Z-2-eq', 'f-Z-3-all',
           'f-Z-3-ax', 'f-Z-3-eq', 'f-chi-0-all', 'f-chi-0-ax', 'f-chi-0-eq', 'f-chi-1-all', 'f-chi-1-ax', 'f-chi-1-eq',
           'f-chi-2-all', 'f-chi-2-ax', 'f-chi-2-eq', 'f-chi-3-all', 'f-chi-3-ax', 'f-chi-3-eq', 'lc-I-0-ax',
           'lc-I-0-eq', 'lc-I-1-ax', 'lc-I-1-eq', 'lc-I-2-ax', 'lc-I-2-eq', 'lc-I-3-ax', 'lc-I-3-eq', 'lc-S-0-ax',
           'lc-S-0-eq', 'lc-S-1-ax', 'lc-S-1-eq', 'lc-S-2-ax', 'lc-S-2-eq', 'lc-S-3-ax', 'lc-S-3-eq', 'lc-T-0-ax',
           'lc-T-0-eq', 'lc-T-1-ax', 'lc-T-1-eq', 'lc-T-2-ax', 'lc-T-2-eq', 'lc-T-3-ax', 'lc-T-3-eq', 'lc-Z-0-ax',
           'lc-Z-0-eq', 'lc-Z-1-ax', 'lc-Z-1-eq', 'lc-Z-2-ax', 'lc-Z-2-eq', 'lc-Z-3-ax', 'lc-Z-3-eq', 'lc-chi-0-ax',
           'lc-chi-0-eq', 'lc-chi-1-ax', 'lc-chi-1-eq', 'lc-chi-2-ax', 'lc-chi-2-eq', 'lc-chi-3-ax', 'lc-chi-3-eq',
           'mc-I-0-all', 'mc-I-1-all', 'mc-I-2-all', 'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all',
           'mc-S-3-all', 'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-Z-0-all', 'mc-Z-1-all',
           'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all', 'mc-chi-2-all', 'mc-chi-3-all']
if tf.__version__ >= tf.__version__ >= '2.0.0':
    tf.compat.v1.disable_eager_execution()


def name_converter(fnames):
    fnames_new = []
    for ii, fname in enumerate(fnames):
        if isRAC(fname):
            fnames_new.append("RACs." + fname)
        elif fname in list(name_converter_dict.keys()):
            fnames_new.append(name_converter_dict[fname])
        else:
            fnames_new.append(fname)
    return fnames_new


def isRAC(fname):
    _fname = fname.split("-")
    if len(_fname) == 4:
        if _fname[0] in ['f', 'lc', 'mc', 'D_lc', 'D_mc'] and _fname[-1] in ['all', 'eq', 'ax']:
            try:
                _ = int(_fname[2])
                return True
            except ValueError:
                return False


def get_vars(predictor):
    fnames = load_ANN_variables(predictor)
    return name_converter(fnames)


def get_label(predictor):
    key = get_key(predictor, suffix="train_y")
    path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.csv')
    _df = pd.read_csv(path_to_file)
    lname = _df.columns.tolist()
    assert len(lname) == 1
    return lname


def get_model_arch(model):
    arch = []
    for layer in model.get_config()['layers']:
        if layer['class_name'] == "Dense":
            arch.append(layer['config']['units'])
    arch = arch[:-1]
    return arch


def get_latest_model(predictor, db):
    cursors = db.models.find({"predictor": predictor})
    model, newest = False, False
    for doc in cursors:
        if model == False:
            model = load_model(str(doc["model"]))
            newest = doc['date']
        else:
            if newest < doc['date']:
                model = load_model(str(doc["model"]))
                newest = doc['date']
    print(('Latest model is trained at: ', newest))
    return model


def extract_data_from_db(predictor, db, collection, constraints,
                         feature_extra=False, target=False):
    print(("Collecting data with constraints: %s..." % constraints))
    df = convert2dataframe(db, collection, constraints=constraints, normalized=True)
    if feature_extra and target:
        print(("Using custom features RACs-180 + ", feature_extra))
        print(("Target property: ", target))
        _fnames = RACs180 + feature_extra
        _fnames = name_converter(_fnames)
        fnames, el = [], []
        for key in _fnames:
            df = df[df[key] != "undef"]
        for f in _fnames:
            std = np.std(df[str(f)].dropna().values)
            if std > 1e-6:
                # print(f, std)
                fnames.append(f)
            else:
                el.append(f)
        print(("features eliminated because of small (<1e-4) std: ", el, len(el)))
        lname = [target]
    elif feature_extra or target:
        raise KeyError("You have to have both or neither <feature_extra> and <target> in input files.")
    else:
        print(("Using predictor-linked features and target.", predictor))
        fnames = get_vars(predictor)
        lname = get_label(predictor)
    print(("features: ", fnames, len(fnames)))
    print(("target: ", lname, len(lname)))
    if lname[0] in list(group_conditions.keys()):
        print(("Paring runs to calculate the new property %s..." % lname))
        df, _ = pairing(df, case=lname[0])
    df_use = df[fnames + lname + ['unique_name']]
    shape = df_use.shape[0]
    df_use = df_use.dropna()
    for key in df_use:
        df_use = df_use[df_use[key] != "undef"]
    print(("data reduce (%d ->  %d) because of NaN." % (shape, df_use.shape[0])))
    return df_use, fnames, lname


def normalize_data(df, fnames, lname, predictor, frac=0.8, name=False):
    np.random.seed(1234)
    X = df[fnames].values
    y = df[lname].values
    if name:
        n = df['unique_name'].values
        X, y, n = sklearn.utils.shuffle(X, y, n)
        n_train, n_test = np.split(n, [int(frac * X.shape[0])])
    else:
        X, y = sklearn.utils.shuffle(X, y)
    X_train, X_test = np.split(X, [int(frac * X.shape[0])])
    y_train, y_test = np.split(y, [int(frac * X.shape[0])])
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_test = x_scaler.transform(X_test)
    if not 'clf' in predictor:
        y_scaler = sklearn.preprocessing.StandardScaler()
        y_scaler.fit(y_train)
        y_train = y_scaler.transform(y_train)
        y_test = y_scaler.transform(y_test)
    print(("mean in target train/test: %f/%f" % (np.mean(y_train), np.mean(y_test))))
    if name:
        return X_train, X_test, y_train, y_test, n_train, n_test
    else:
        return X_train, X_test, y_train, y_test


def train_model(predictor, db,
                X_train, X_test, y_train, y_test,
                epochs=1000, batch_size=32, hyperopt=False,
                initialize_weight=True, hyperopt_step=100,
                load_latest_model=False, fix_architecture=False):
    regression = False if 'clf' in predictor else True
    if not hyperopt:  # load molSimplify model and training parameters directly.
        model = load_keras_ann(predictor)
        if initialize_weight:
            print("Initializing weights for the final model training...")
            model = initialize_model_weights(model)
        best_params = {"epochs": epochs, "batch_size": batch_size}
    else:
        if fix_architecture:  # HyperOpt, with a fix architecture but flexible training parameters.
            if not load_latest_model:  # Use molSimplify pretrained model architecture.
                model = load_keras_ann(predictor)
            else:  # Use the lastest model achitecture.
                model = get_latest_model(predictor, db)
            if not model == False:
                if initialize_weight:
                    print("Initializing weights for the final model training...")
                    model = initialize_model_weights(model)
                arch = get_model_arch(model)
                print(('HyperOpt with a FIXED architecture: ', arch))
                best_params = optimize(X=X_train, y=y_train,
                                       regression=regression,
                                       hyperopt_step=hyperopt_step,
                                       arch=arch,
                                       epochs=epochs)
            else:
                print(("No existing model found in db.models matching the predictor: ", predictor))
                print('HyperOpt EVERYTHING...')
                best_params = optimize(X=X_train, y=y_train,
                                       regression=regression,
                                       hyperopt_step=hyperopt_step,
                                       epochs=epochs)
        else:  # HyperOpt both model architecture and training parameters.
            print('HyperOpt EVERYTHING...')
            best_params = optimize(X=X_train, y=y_train,
                                   regression=regression,
                                   hyperopt_step=hyperopt_step,
                                   epochs=epochs)
        print(("best hyperparams: ", best_params))
        model = build_ANN(best_params, input_len=X_train.shape[-1], regression=regression)
        batch_size = best_params['batch_size']
        epochs = best_params['epochs']
        print(("epochs: %d, batch_size: %d" % (epochs, batch_size)))
    history = model.fit(X_train, y_train, epochs=epochs, verbose=2, batch_size=batch_size)
    results = model.evaluate(X_train, y_train)
    res_dict_train = {}
    for ii, key in enumerate(model.metrics_names):
        res_dict_train.update({key: results[ii]})
    results = model.evaluate(X_test, y_test)
    res_dict_test = {}
    for ii, key in enumerate(model.metrics_names):
        res_dict_test.update({key: results[ii]})
    if not regression:
        res_dict_train.update({"auc": cal_auc(model, X_train, y_train)})
        res_dict_test.update({"auc": cal_auc(model, X_test, y_test)})
    print(("train result: ", res_dict_train))
    print(("test reulst: ", res_dict_test))
    return model, history, res_dict_train, res_dict_test, best_params


def retrain(predictor, user, pwd,
            database, collection, collection_model,
            host="localhost", port=27017, auth=True,
            constraints=False, frac=0.8, epochs=1000,
            batch_size=32, force_push=False,
            hyperopt=False, tag=False,
            feature_extra=False, target=False,
            initialize_weight=True, hyperopt_step=100,
            load_latest_model=False, fix_architecture=False):
    db = connect2db(user, pwd, host, port, database, auth)
    dbquery_time = datetime.now()
    df, fnames, lname = extract_data_from_db(predictor, db, collection,
                                             constraints=constraints,
                                             feature_extra=feature_extra,
                                             target=target)
    X_train, X_test, y_train, y_test, n_train, n_test = normalize_data(df, fnames, lname, predictor,
                                                                       frac=frac, name=True)
    model, history, res_dict_train, res_dict_test, best_params = train_model(predictor, db,
                                                                             X_train, X_test,
                                                                             y_train, y_test,
                                                                             epochs=epochs,
                                                                             batch_size=batch_size,
                                                                             hyperopt=hyperopt,
                                                                             initialize_weight=initialize_weight,
                                                                             hyperopt_step=hyperopt_step,
                                                                             load_latest_model=load_latest_model,
                                                                             fix_architecture=fix_architecture)
    model_dict = {}
    model_dict.update({"predictor": predictor,
                       "constraints": str(constraints),
                       "hyperopt": hyperopt,
                       "history": {k: [float(ele) for ele in history.history[k]] for k in history.history},
                       "hyperparams": best_params,
                       "dbquery_time": dbquery_time,
                       "score_train": {k: float(res_dict_train[k]) for k in res_dict_train},
                       "score_test": {k: float(res_dict_test[k]) for k in res_dict_test},
                       "name_train": n_train.tolist(),
                       "name_test": n_test.tolist(),
                       "target_train": y_train.tolist(),
                       "target_test": y_test.tolist(),
                       "pred_train": model.predict(X_train).tolist(),
                       "pred_test": model.predict(X_test).tolist(),
                       "len_train": y_train.shape[0],
                       "len_test": y_test.shape[0],
                       "len_tot": y_train.shape[0] + y_test.shape[0],
                       "features": fnames,
                       "target": lname,
                       "force_push": force_push,
                       "tag": tag,
                       "initialize_weight": initialize_weight,
                       "hyperopt_step": hyperopt_step,
                       "load_latest_model": load_latest_model,
                       "fix_architecture": fix_architecture,
                       })
    push_models(model, model_dict,
                database, collection_model,
                user=user, pwd=pwd,
                host=host, port=port,
                auth=auth)


def check_retrain_inputs(args_dict):
    keys_required = ["predictor", "database", "collection", "collection_model", "user", "pwd"]
    if not set(keys_required) <= set(args_dict.keys()):
        raise KeyError("missing necessary keys to retrain a model. Required inputs are: ", keys_required)
    return args_dict


def retrain_and_push(args_dict):
    args_dict = check_retrain_inputs(args_dict)
    default_args = {"host": "localhost", "port": 27017, "auth": True,
                    "constraints": False, "frac": 0.8, "epochs": 1000,
                    "batch_size": 32, "force_push": False,
                    "hyperopt": False, "tag": False,
                    "feature_extra": False, "target": False,
                    "initialize_weight": True, "hyperopt_step": 100,
                    "load_latest_model": False, "fix_architecture": False}
    for key in args_dict:
        print((key, args_dict[key]))
        globals().update({key: args_dict[key]})
    for key in default_args:
        if not key in list(args_dict.keys()):
            globals().update({key: default_args[key]})
    retrain(predictor, user, pwd,
            database, collection, collection_model,
            host, port, auth,
            constraints, frac, epochs,
            batch_size, force_push,
            hyperopt, tag,
            feature_extra, target,
            initialize_weight, hyperopt_step,
            load_latest_model, fix_architecture)
