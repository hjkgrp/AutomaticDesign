import os
from pymongo import MongoClient
import pandas as pd
import pickle
from molSimplifyAD.dbclass_mongo import tmcMongo
from molSimplifyAD.ga_check_jobs import check_all_current_convergence
from molSimplifyAD.ga_tools import isKeyword


def check_repeated(db, collection, tmc):
    doc = db[collection].find_one(tmc.id_doc)
    if doc == None:
        repeated = False
    else:
        repeated = True
    return repeated, doc


def query_to_db(db, collection, constraints):
    return db[collection].find(constraints)


def query_one(db, collection, constraints):
    return db[collection].find_one(constraints)


def insert(db, collection, tmc):
    repeated, _tmc = check_repeated(db, collection, tmc)
    inserted = False
    if not repeated:
        db[collection].insert_one(tmc.document)
        inserted = True
    else:
        print("existed: ", _tmc.id_doc)
    return inserted


def connect_to_db(user, pwd, host, port, database):
    cstring = "mongodb://%s:%s@%s:%s/%s/" % (user, pwd, host, port, database)
    client = MongoClient(cstring)
    db = client[database]
    return db


def deserialize_dftrun_from_db(tmc):
    if tmc["dftrun"]:
        this_run = pickle.loads(tmc["dftrun"])
    else:
        this_run = False
    return this_run


def convert_to_dataframe(db, collection, constraints=False):
    if constraints:
        cursor = query_to_db(db, collection, constraints)
    else:
        cursor = db[collection]
    df = pd.DataFrame(list(cursor))
    return df


def push_to_db(user, pwd, host, port, database, collection):
    if not os.path.isfile(".madconfig") or not isKeyword('tag'):
        raise ValueError("This is not a mAD folder with a tag to push.")
    tag = str(isKeyword('tag'))
    db = connect_to_db(user, pwd, host, port, database)
    colls = db.list_collection_names()
    if not collection in colls:
        finish = False
        while not finish:
            print("Collection %s is not currently in this database. Are you sure you want to create a new collection? (y/n)")
            _in = str(input())
            if _in == "y":
                finish = True
            elif _in == "n":
                print("Quit. Have a nive day.")
                quit()
            else:
                finish = False
    print('Warning: db push is enabled, attempting commit with tag ' + str(isKeyword('tag')))
    _, all_runs = check_all_current_convergence()
    count = 0
    for this_run in all_runs.values():
        this_tmc = tmcMongo(this_run=this_run, tag=tag)
        insetred = insert(db, collection, this_tmc)
        if insetred:
            count += 1
    print("add %d entries in the %s['%s']." % (count, database, collection))
