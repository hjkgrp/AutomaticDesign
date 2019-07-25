import numpy as np
import pandas as pd
import sklearn.preprocessing
import sklearn.utils
from molSimplifyAD.utils.pymongo_tools import convert2dataframe, connect2db, push_models
from pkg_resources import resource_filename, Requirement
from molSimplify.python_nn.tf_ANN import get_key, load_ANN_variables, load_keras_ann, initialize_model_weights

name_converter_dict = {"oxstate": "ox", "spinmult": "spin", "charge_lig": "ligcharge"}


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


def extract_data_from_db(predictor, db, collection, constraints):
    print(("Collecting data with constraints: %s..." % constraints))
    df = convert2dataframe(db, collection, constraints=constraints, normalized=True)
    fnames = get_vars(predictor)
    lname = get_label(predictor)
    df_use = df[fnames + lname]
    shape = df_use.shape[0]
    df_use = df_use.dropna()
    for key in df_use:
        df_use = df_use[df_use[key] != "undef"]
    print("data reduce (%d ->  %d) because of NaN." % (shape, df_use.shape[0]))
    return df_use, fnames, lname


def normalize_data(df, fnames, lname, predictor, frac=0.8):
    np.random.seed(1234)
    X = df[fnames].values
    y = df[lname].values
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
    return X_train, X_test, y_train, y_test


def train_model(predictor, X_train, X_test, y_train, y_test, epochs=1000, batch_size=32):
    model = load_keras_ann(predictor)
    print("Initializing weights...")
    model = initialize_model_weights(model)
    history = model.fit(X_train, y_train, epochs=epochs, verbose=1, batch_size=batch_size)
    loss, metrics = model.evaluate(X_test, y_test)
    print(("loss: ", loss))
    if not 'clf' in predictor:
        print(("mae: ", metrics))
    else:
        print("accuracy: ", metrics)
    return model, history


def retrain(predictor, user, pwd,
            database, collection, collection_model,
            host="localhost", port=27017, auth=True,
            constraints=False, frac=0.8, epochs=1000,
            batch_size=32, force_push=False):
    db = connect2db(user, pwd, host, port, database, auth)
    df, fnames, lname = extract_data_from_db(predictor, db, collection, constraints=constraints)
    X_train, X_test, y_train, y_test = normalize_data(df, fnames, lname, predictor, frac=frac)
    model, history = train_model(predictor, X_train, X_test, y_train, y_test,
                                 epochs=epochs, batch_size=batch_size)
    model_dict = {}
    model_dict.update({"predictor": predictor})
    model_dict.update({"constraints": str(constraints)})
    model_dict.update({"history": history.history})
    model_dict.update({"hyperparams": {"epochs": epochs, "batch_size": batch_size}})
    model_dict.update({"score_train": model.evaluate(X_train, y_train)[-1],
                       "score_test": model.evaluate(X_test, y_test)[-1],
                       "target_train": y_train.tolist(),
                       "target_test": y_test.tolist(),
                       "pred_train": model.predict(X_train).tolist(),
                       "pred_test": model.predict(X_test).tolist(),
                       "len_train": y_train.shape[0],
                       "len_test": y_test.shape[0],
                       "len_tot": y_train.shape[0] + y_test.shape[0],
                       "force_push": force_push
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
                    "batch_size": 32, "force_push": False}
    for key in args_dict:
        print(key, args_dict[key])
        globals().update({key: args_dict[key]})
    for key in default_args:
        if not key in args_dict.keys():
            globals().update({key: default_args[key]})
    retrain(predictor, user, pwd,
            database, collection, collection_model,
            host, port, auth,
            constraints, frac, epochs,
            batch_size, force_push)
