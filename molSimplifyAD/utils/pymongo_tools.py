import os
import sys
from pymongo import MongoClient
import pandas as pd
from pandas.io.json import json_normalize
import pickle
from molSimplifyAD.dbclass_mongo import tmcMongo, mongo_attr_id
from molSimplifyAD.ga_check_jobs import check_all_current_convergence
from molSimplifyAD.ga_tools import isKeyword


def check_repeated(db, collection, tmc):
    doc = db[collection].find_one(tmc.id_doc)
    if doc == None:
        repeated = False
    else:
        repeated = True
    return repeated, doc


def query_db(db, collection, constraints):
    return db[collection].find(constraints)


def query_one(db, collection, constraints):
    return db[collection].find_one(constraints)


def insert(db, collection, tmc):
    repeated, _tmcdoc = check_repeated(db, collection, tmc)
    inserted = False
    if not repeated:
        db[collection].insert_one(tmc.document)
        inserted = True
    else:
        this_tmc = tmcMongo(document=_tmcdoc, tag="undef", subtag='undef')
        print("existed: ", this_tmc.id_doc)
    return inserted


def convert2dataframe(db, collection, constraints=False, dropcols=["dftrun"], directload=False, normalized=False):
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


def push2db(database, collection, user=False, pwd=False, host="localhost", port=27017, auth=False,
            all_runs_pickle=False, tag=False, subtag=False):
    if not tag:
        if not os.path.isfile(".madconfig") or not isKeyword('tag'):
            raise ValueError("This is not a mAD folder with a tag to push.")
        tag = str(isKeyword('tag'))
    if not subtag:
        if not os.path.isfile(".madconfig") or not isKeyword('subtag'):
            raise ValueError("This is not a mAD folder with a subtag to push.")
        subtag = str(isKeyword('subtag'))
    db = connect2db(user, pwd, host, port, database, auth)
    colls = db.list_collection_names()
    if not collection in colls:
        finish = False
        while not finish:
            print(
                "Collection %s is not currently in this database. Are you sure you want to create a new collection? (y/n)" % collection)
            _in = raw_input()
            if _in == "y":
                finish = True
            elif _in == "n":
                print("Quit. Have a nive day.")
                quit()
            else:
                finish = False
    print('Warning: db push is enabled, attempting commit with tag: %s, subtag: %s' % (tag, subtag))
    if not all_runs_pickle:
        _, all_runs = check_all_current_convergence()
    else:
        all_runs = pickle.load(open(all_runs_pickle, "rb"))
        print("DFTruns loaded from %s." % all_runs_pickle)
    count = 0
    for this_run in all_runs.values():
        print("complex: ", this_run.name)
        if sys.getsizeof(pickle.dumps(this_run)) * 1. / 10 ** 6 > 16.7:
            print(
                "DFTrun too large. Deleting wavefunction binary. Only the path of wavefunction files is hold by DFTrun.")
            _this_tmc = tmcMongo(this_run=this_run, tag=tag, subtag=subtag)
            this_run = _this_tmc.write_wfn(this_run)
        this_tmc = tmcMongo(this_run=this_run, tag=tag, subtag=subtag)
        insetred = insert(db, collection, this_tmc)
        if insetred:
            count += 1
    print("add %d entries in the %s['%s']." % (count, database, collection))


def unique_name(tmcdoc):
    this_tmc = tmcMongo(document=tmcdoc, tag=tmcdoc['tag'], subtag=tmcdoc['subtag'])
    name_ele = []
    for key in mongo_attr_id:
        name_ele.append(key)
        name_ele.append(str(this_tmc.id_doc[key]))
    name = '_'.join(name_ele)
    return name
