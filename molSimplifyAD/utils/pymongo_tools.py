import os
import sys
import pickle
import time
import pandas as pd
import pymongo
from pandas.io.json import json_normalize
from pymongo import MongoClient
from molSimplifyAD.dbclass_mongo import tmcMongo, tmcActLearn, mongo_attr_id, mongo_not_web
from molSimplifyAD.mlclass_mongo import modelActLearn, modelMongo
from molSimplifyAD.ga_tools import isKeyword


def check_repeated(db, collection, tmc):
    doc = db[collection].find_one(tmc.id_doc)
    if doc == None:
        repeated = False
    else:
        repeated = True
    return repeated, doc


def query_db(db, collection, constraints):
    '''
    Query the databse.

    :param db: mongo database instance.
    :param collection: name of a collection.
    :param constraints: a dictionary of conditions for query.
    :return: a cursor.
    '''
    return db[collection].find(constraints)


def query_one(db, collection, constraints):
    return db[collection].find_one(constraints)


def query_lowestE_converged(db, collection, constraints):
    cursor = query_db(db, collection, constraints)
    minE = 100000
    tmcdoc = None
    for _tmcdoc in cursor:
        if _tmcdoc["converged"] and _tmcdoc["energy"] < minE:
            minE = _tmcdoc["energy"]
            tmcdoc = _tmcdoc
    return tmcdoc


def insert(db, collection, tmc):
    repeated, _tmcdoc = check_repeated(db, collection, tmc)
    inserted = False
    if not repeated:
        db[collection].insert_one(tmc.document)
        inserted = True
    else:
        print("existed: ", _tmcdoc["unique_name"])
        print("merging....")
        _tmc = tmcMongo(document=_tmcdoc, tag=_tmcdoc["tag"], subtag=_tmcdoc["subtag"],
                        publication=_tmcdoc["publication"])
        merge_documents(db, collection,
                        doc1=_tmcdoc, doc2=tmc.document,
                        update_fields=tmc.update_fields)
        merge_dftruns(dftrun1=_tmc.this_run,
                      dftrun2=tmc.this_run,
                      update_fields=tmc.update_fields)
        _tmc.write_dftrun()
    return inserted


def merge_documents(db, collection, doc1, doc2, update_fields):
    for key in doc2:
        if not key in doc1 or key in update_fields:
            doc1.update({key: doc2[key]})
    db[collection].replace_one({"_id": doc1["_id"]}, doc1)


def merge_dftruns(dftrun1, dftrun2, update_fields):
    for attr, val in dftrun2.__dict__.items():
        if not attr in dftrun1.__dict__ or attr in update_fields:
            setattr(dftrun1, attr, val)


def convert2dataframe(db, collection,
                      constraints=False, dropcols=["dftrun"],
                      directload=False, normalized=False):
    '''
    Converts a collection into pandas dataframe.

    :param db: mongo database instance.
    :param collection: name of a collection.
    :param constraints: a dictionary of conditions for query.
    :param dropcols: a list of columns to drop when converting a collection to pandas df.
    :param directload: whether to load the whole collection at one time. Might be slow if a collection is large.
    :param normalized: whether to json normalize each column of the collection.
    :return: a pandas dataframe
    '''
    if constraints:
        cursor = query_db(db, collection, constraints)
    else:
        cursor = db[collection].find()
    if directload:
        if not normalized:
            df = pd.DataFrame(list(cursor))
        else:
            df = json_normalize(list(cursor))
    else:
        df = iterator2dataframes_withoutruns(cursor, 100, dropcols, normalized)
    return df


def iterator2dataframes_withoutruns(iterator, chunk_size, dropcols=["dftrun"], normalized=False):
    records = []
    frames = []
    for i, record in enumerate(iterator):
        records.append(record)
        if i % chunk_size == chunk_size - 1:
            if not normalized:
                frames.append(pd.DataFrame(records).drop(dropcols, axis=1))
            else:
                frames.append(json_normalize(records).drop(dropcols, axis=1))
            records = []
    if records:
        if not normalized:
            frames.append(pd.DataFrame(records).drop(dropcols, axis=1))
        else:
            frames.append(json_normalize(records).drop(dropcols, axis=1))
    return pd.concat(frames)


def deserialize_dftrun_from_db(tmc):
    if tmc["dftrun"]:
        this_run = pickle.loads(tmc["dftrun"])
    else:
        this_run = False
    return this_run


def connect2db(user, pwd, host, port, database, auth):
    if not auth:
        client = MongoClient()
    else:
        cstr = "mongodb://%s:%s@%s:%d/%s" % (user, pwd, host, port, database)
        client = MongoClient(cstr)
    db = client[database]
    return db


def ensure_collection(db, collection):
    colls = db.list_collection_names()
    if not collection in colls:
        finish = False
        while not finish:
            print("Collection %s does not exist. Create a new collection? (y/n)" % collection)
            _in = raw_input()
            if _in == "y":
                finish = True
            elif _in == "n":
                print("Quit. Have a nive day.")
                quit()
            else:
                finish = False


def push2db(database, collection, tag, subtag, publication=False,
            user=False, pwd=False, host="localhost", port=27017, auth=False,
            all_runs_pickle=False, update_fields=False):
    '''
    Push data to MongoDB.

    :param db: mongo database instance.
    :param collection: name of a collection.
    :param user: username.
    :param pwd: password.
    :param host: IP address of the MongoDB.
    :param port: port to connect. 
    :param auth: whether authentication is required to connect to the MongoDB.
    :param all_runs_pickle: whether to push from a pickle file of a list of DFTrun objects. If False, will run
    check_all_current_convergence() in your mAD folder.
    :param tag: tag of your data. Recommend to use the project name (may related to the name of your paper).
    :param subtag: Recommend as the name of your mAD folder (considering we may have many mAD folders for each project).

    :return: None
    '''
    from molSimplifyAD.ga_check_jobs import check_all_current_convergence
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    print('db push is enabled, attempting commit with tag: %s, subtag: %s to %s' % (tag, subtag, collection))
    if not all_runs_pickle:
        _, all_runs = check_all_current_convergence(post_all=True)
    else:
        all_runs = pickle.load(open(all_runs_pickle, "rb"))
        print("DFTruns loaded from %s." % all_runs_pickle)
    print("number of DFTruns exists: ", len(all_runs))
    count = 0
    merged = 0
    for this_run in all_runs.values():
        print("adding complex: ", this_run.name)
        print("converged:", this_run.converged, "geo_flag:", this_run.geo_flag,
              "ss_flag: ", this_run.ss_flag, "metal_spin_flag: ", this_run.metal_spin_flag)
        this_tmc = tmcMongo(this_run=this_run, tag=tag, subtag=subtag,
                            publication=publication, update_fields=update_fields)
        _s = time.time()
        insetred = insert(db, collection, this_tmc)
        print("elapse: ", time.time() - _s)
        if insetred:
            count += 1
        else:
            merged += 1
    print("add %d entries in the %s['%s']." % (count, database, collection))
    print("merge %d entries in the %s['%s']." % (merged, database, collection))
    print("creating index...")
    db[collection].create_index([("metal", pymongo.ASCENDING),
                                 ("ox", pymongo.ASCENDING),
                                 ("spin", pymongo.ASCENDING),
                                 ("converged", pymongo.ASCENDING),
                                 ("alpha", pymongo.ASCENDING),
                                 ("lig1", pymongo.ASCENDING),
                                 ("lig5", pymongo.ASCENDING),
                                 ("lig6", pymongo.ASCENDING)
                                 ])


def push_complex_actlearn(step, all_complexes, database, collection,
                          user=False, pwd=False, host="localhost", port=27017,
                          auth=False):
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    print('db push is enabled, attempting commit to' % collection)
    print("number of complexes to push: ", len(all_complexes))
    count = 0
    merged = 0
    for this_complex in all_complexes:
        print("adding complex: ", this_complex["dftrun"].name)
        this_tmc = tmcActLearn(step=step,
                               is_training=this_complex["is_training"],
                               status_flag=this_complex["status_flag"],
                               target=this_complex["target"],
                               descriptors=this_complex["descriptors"],
                               this_run=this_complex["dftrun"]
                               )
        _s = time.time()
        repeated, _tmcdoc = check_repeated(db, collection, this_tmc)
        if not repeated:
            db[collection].insert_one(this_tmc.document)
            count += 1
        else:
            merged += 1
        print("elapse: ", time.time() - _s)
    print("add %d entries in the %s['%s']." % (count, database, collection))
    print("merge %d entries in the %s['%s']." % (merged, database, collection))
    print("creating index...")
    db[collection].create_index([("step", pymongo.ASCENDING),
                                 ("is_training", pymongo.ASCENDING),
                                 ("status_flag", pymongo.ASCENDING),
                                 ("ox", pymongo.ASCENDING),
                                 ("spin", pymongo.ASCENDING),
                                 ("converged", pymongo.ASCENDING),
                                 ("alpha", pymongo.ASCENDING),
                                 ("lig1", pymongo.ASCENDING),
                                 ("lig5", pymongo.ASCENDING),
                                 ("lig6", pymongo.ASCENDING)
                                 ])
    if not merged == 0:
        print("=====WARNING====")
        print("Duplicate complexes(%d) occure in the active learning mode. Should never happen." % merged)


def push_models(model, model_dict, database, collection,
                user=False, pwd=False,
                host="localhost", port=27017,
                auth=False):
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    this_model = modelMongo(model=model, **model_dict)
    if not query_one(db, collection,
                     constraints={"predictor": this_model.predictor, "len_tot": this_model.len_tot}) == None:
        print("A model of step %d has already existed.")
    elif this_model.force_push:
        print("force_push is truned on. pushing...")
        db[collection].insert_one(this_model.document)
    else:
        db[collection].insert_one(this_model.document)


def push_moldels_actlearn(step, model, database, collection,
                          user=False, pwd=False,
                          host="localhost", port=27017,
                          auth=False):
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    actlearn_model = modelActLearn(step=step, model=model)
    if not query_one(db, collection, constraints={"step": step}) == None:
        print("A model of step %d has already existed.")
    else:
        db[collection].insert_one(actlearn_model.document)
