import os
import sys
import pickle
import subprocess
import time
import datetime
import pandas as pd
import numpy as np
import pymongo
import itertools
from pandas import json_normalize
from pymongo import MongoClient
from molSimplifyAD.dbclass_mongo import tmcMongo, tmcActLearn, mongo_attr_id, mongo_not_web, SP_keys
from molSimplifyAD.mlclass_mongo import modelActLearn, modelMongo, modelPublished
from molSimplifyAD.ga_tools import isKeyword


def check_repeated(db, collection, tmc):
    doc = db[collection].find_one(tmc.id_doc)
    if doc == None:
        repeated = False
    else:
        repeated = True
    return repeated, doc


def dump_databse(database_name="tmc", outpath='/userarchive/db_backup',
                 user=False, pwd=False):
    now = datetime.datetime.now()
    run_cmd = 'mongodump -d %s -o %s -u %s -p %s' % (database_name, outpath, user, pwd)
    q = subprocess.Popen(run_cmd, shell=True, stdout=subprocess.PIPE)
    print("Dumping database....")
    ll = q.communicate()[0].decode("utf-8")
    print(ll)
    print("Done.")
    with open(outpath + '/dump.log', "a") as fo:
        fo.write("dumping database %s to path %s at time %s.\n" % (database_name, outpath, str(now)))


def query_db_easy(user, pwd, collection='oct', constraints={}, max_entries=None, dropcols=[], host='localhost',
                  port=27017, database='tmc', auth=True, loud=True):
    '''
    Query the database for `max_entries` entries using constraints in `constraints`. Queries the octahedral transition metal complex database by default.

    Example usage: results_dataframe = query_db_easy(username, password, constraints={"geo_flag": 1}, max_entries=100, dropcols=['dftrun'], loud=False)
    Returns: Pandas dataframe with database query results.
    Example constraints string: constraints = {"metal": {"$in": ['fe', 'co']}, "tot_time": {"$gt": 0}, "ox": 2}
    '''

    start_time = time.time()
    db = connect2db(user=user, pwd=pwd, host=host, port=port, database=database, auth=auth)

    if collection == None:
        raise ValueError(
            'Argument `collection` not specified; you can use the following collections: %s' % db.list_collection_names())
    if loud:
        if constraints == {}:
            print("Warning: argument `constraints` is {} by defaultl. No constraints are being applied.")
        if max_entries == None:
            print(
                "Warning: argument `max_entries` is None by default, so no limit has been set on the number of entries to pull.")
        if dropcols == []:
            print("Warning: argument `dropcols` is [] by default, so no columns will be dropped.")

    est_num_docs = db.oct.count_documents(constraints)
    num_docs = est_num_docs if max_entries == None else min(est_num_docs, max_entries)
    if loud:
        print("Estimated number of entries to be pulled: %s" % num_docs)
        print("Estimated pull time: %s sec" % (2.36e-4 * num_docs))  # Empirical fit, as of 2020-01-16
        print("Estimated pull size: %s MB" % (3.53e-2 * num_docs))  # Empirical fit, as of 2020-01-16

    if dropcols == []:
        cursor = db[collection].find(constraints)
    else:
        excluded_fields = {dropcol: 0 for dropcol in dropcols}
        cursor = db[collection].find(constraints, excluded_fields)
    cursor_iterator = itertools.islice(cursor, max_entries)
    results = pd.DataFrame(cursor_iterator)
    if loud:
        print("Pull finished. Number of entries: %s" % len(results))
        print("Time for pull: %s sec" % (time.time() - start_time))
    return results


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


def count_find(cursor):
    count = 0
    for _ in cursor:
        count += 1
    return count


def query_lowestE_converged(db, collection, constraints):
    cursor = query_db(db, collection, constraints)
    minE = 100000
    tmcdoc = None
    for _tmcdoc in cursor:
        if _tmcdoc["converged"] and _tmcdoc["energy"] < minE:
            minE = _tmcdoc["energy"]
            tmcdoc = _tmcdoc
    return tmcdoc


def merge_document_with_new(db, doc1, doc_new, update_fields):
    for k in doc1:
        if not k in update_fields:
            doc_new.update({k: doc1[k]})
    for k in doc_new:
        if k in doc1 and not k in update_fields:
            if not str(doc1[k]) == str(doc_new[k]):
                raise ValueError("Error on %s. If you want to overwrite existing key-value pairs in the document, add it to the <update_fields>"%k)
    for k in doc1:
        if not k in doc_new:
            raise ValueError("Error on %s. You are removing an existing key-value pair in the document. Please think twice before you do this."%k)
    print("replacing doc without changing dftrun: ", doc1['_id'], doc1['unique_name'])
    db.oct.find_one_and_replace({"_id": doc1['_id']}, doc_new)


def insert(db, collection, tmc, debug=True):
    repeated, _tmcdoc = check_repeated(db, collection, tmc)
    inserted = False
    if not repeated:
        db[collection].insert_one(tmc.document)
        inserted = True
    else:
        if debug:
            print(("existed: ", _tmcdoc["unique_name"]))
            print("merging....")
        try:
            _tmc = tmcMongo(document=_tmcdoc, tag=_tmcdoc["tag"], subtag=_tmcdoc["subtag"],
                            publication=_tmcdoc["publication"])
            converted = True
        except ValueError:
            converted = False
            merge_document_with_new(db, doc1=_tmcdoc, doc_new=tmc.document, update_fields=tmc.update_fields)
        if converted:
            merge_documents(db, collection,
                            doc1=_tmcdoc, doc2=tmc.document,
                            update_fields=tmc.update_fields)
            merge_dftruns(dftrun1=_tmc.this_run,
                          dftrun2=tmc.this_run,
                          update_fields=tmc.update_fields)
            _tmc.write_dftrun(force=True)
    return inserted


def merge_documents(db, collection, doc1, doc2, update_fields):
    for key in doc2:
        if not key in doc1:
            doc1.update({key: doc2[key]})
        elif key in update_fields:
            if key in SP_keys and type(doc1[key]) == dict and type(doc2[key]) == dict:
                doc1[key].update(doc2[key])
            else:
                doc1.update({key: doc2[key]})
    db[collection].replace_one({"_id": doc1["_id"]}, doc1)


def merge_dftruns(dftrun1, dftrun2, update_fields):
    for attr, val in list(dftrun2.__dict__.items()):
        if not attr in dftrun1.__dict__:
            setattr(dftrun1, attr, val)
        elif attr in update_fields:
            if attr in SP_keys and type(getattr(dftrun1, attr)) == dict and type(getattr(dftrun2, attr)) == dict:
                getattr(dftrun1, attr).update(getattr(dftrun2, attr))
            else:
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
        df = iterator2dataframes_withoutruns(cursor, 256, dropcols, normalized)
    return df


def convert2readablecsv(db, collection,
                        constraints=False, dropRACs=False,
                        directload=False, normalized=False,
                        ):
    dropcols = ['status', 'grad_max_hist', 'tot_step', 'dipole_vec', 'displace_rms_hist',
                'functional', 'liglist', 'date', 'init_ligand_symmetry', 'outpath',
                'init_geo', 'old_dynamic_feature', 'initRACs', 'prog_geo', 'dftrun',
                'geo_check_metrics', 'author', 'wavefunction', 'functionals', 'molden',
                'terachem_version', 'step_qual_hist', 'basis', 'tot_time', 'solvent',
                'e_delta_hist', 'scrpath', 'init_mol_graph_det', 'geo_check_metrics_prog',
                'trust_radius_hist', 'dynamic_feature', 'expected_delE_hist', 'e_hist',
                'name', 'geo_check_dict', 'geo_opt', 'grad_rms_hist', 'opt_geo', 'd3opt_flag',
                'csd_doi', 'csd_mol2string', 'dupe_refcode_plus_list', 'diople_vec',
                'is_csd_init_geo', 'iscsd', 'isMutation', 'refcode',
                ]
    if dropRACs:
        dropcols += ["RACs", "lacRACs"]
    df = convert2dataframe(db, collection,
                           constraints=constraints, dropcols=dropcols,
                           directload=directload, normalized=False)
    ordered_cols = ["chemical_name", "metal", "ox", "spin", "ligstr", "charge", "ligcharge", "alpha",
                    'refcode_plus', "tag", "subtag", "doi", 'publication', 'energy']
    if normalized:
        df = json_normalize(df.to_dict("records"))
    cols = list(df.columns)
    cols = ordered_cols + sorted(list(set(cols) - set(ordered_cols)))
    df = df[cols]
    return df


def iterator2dataframes_withoutruns(iterator, chunk_size, dropcols=["dftrun"], normalized=False):
    records = []
    frames = []
    for i, record in enumerate(iterator):
        records.append(record)
        if i % chunk_size == chunk_size - 1:
            if not normalized:
                frames.append(pd.DataFrame(records).drop(dropcols, axis=1, errors='ignore'))
            else:
                frames.append(json_normalize(records).drop(dropcols, axis=1, errors='ignore'))
            records = []
    if records:
        if not normalized:
            frames.append(pd.DataFrame(records).drop(dropcols, axis=1, errors='ignore'))
        else:
            frames.append(json_normalize(records).drop(dropcols, axis=1, errors='ignore'))
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
            print(("Collection %s does not exist. Create a new collection? (y/n)" % collection))
            if sys.version_info[0] < 3:
                _in = raw_input()
            else:
                _in = input()
            if _in == "y":
                finish = True
            elif _in == "n":
                print("Quit. Have a nice day.")
                quit()
            else:
                finish = False


def push2db(database, collection, tag, subtag, publication=False,
            user=False, pwd=False, host="localhost", port=27017, auth=False,
            all_runs_pickle=False, update_fields=False, all_runs_list=False,
            outpath='/userarchive/db_backup'):
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
    print(('db push is enabled, attempting commit with tag: %s, subtag: %s to %s' % (tag, subtag, collection)))
    if not all_runs_pickle:
        if not all_runs_list == False:
            all_runs = all_runs_list
        else:
            _, all_runs, _ = check_all_current_convergence(post_all=True)
    else:
        all_runs = pickle.load(open(all_runs_pickle, "rb"))
        print(("DFTruns loaded from %s." % all_runs_pickle))
    print(("number of DFTruns exists: ", len(all_runs)))
    count = 0
    merged = 0
    for this_run in list(all_runs.values()):
        print(("adding complex: ", this_run.name))
        print(("converged:", this_run.converged, "geo_flag:", this_run.geo_flag,
               "ss_flag: ", this_run.ss_flag, "metal_spin_flag: ", this_run.metal_spin_flag))
        this_tmc = tmcMongo(this_run=this_run, tag=tag, subtag=subtag,
                            publication=publication, update_fields=update_fields)
        _s = time.time()
        insetred = insert(db, collection, this_tmc)
        print(("elapse: ", time.time() - _s))
        if insetred:
            count += 1
        else:
            merged += 1
    print(("add %d entries in the %s['%s']." % (count, database, collection)))
    print(("merge %d entries in the %s['%s']." % (merged, database, collection)))
    print("creating index...")
    db[collection].create_index([("metal", pymongo.ASCENDING),
                                 ("ox", pymongo.ASCENDING),
                                 ("spin", pymongo.ASCENDING),
                                 ("converged", pymongo.ASCENDING),
                                 ("alpha", pymongo.ASCENDING),
                                 ("lig1", pymongo.ASCENDING),
                                 ("lig5", pymongo.ASCENDING),
                                 ("lig6", pymongo.ASCENDING),
                                 ("status", pymongo.ASCENDING),
                                 ("geo_flag", pymongo.ASCENDING),
                                 ("ss_flag", pymongo.ASCENDING)
                                 ])
    dump_databse(database_name=database,
                 outpath=outpath,
                 user=user, pwd=pwd)


def push_complex_actlearn(step, all_complexes, database, collection,
                          user=False, pwd=False, host="localhost", port=27017,
                          auth=False, update_fields=False,
                          outpath='/userarchive/db_backup'):
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    print(('db push is enabled, attempting commit to ', collection))
    print(("number of complexes to push: ", len(all_complexes)))
    count = 0
    merged = 0
    for this_complex in all_complexes:
        print(("adding complex: ", this_complex["dftrun"].name))
        this_tmc = tmcActLearn(step=step,
                               is_training=this_complex["is_training"],
                               status_flag=this_complex["status_flag"],
                               target=this_complex["target"],
                               descriptors=this_complex["descriptors"],
                               this_run=this_complex["dftrun"],
                               update_fields=update_fields
                               )
        _s = time.time()
        repeated, _tmcdoc = check_repeated(db, collection, this_tmc)
        if not repeated:
            db[collection].insert_one(this_tmc.document)
            count += 1
        else:
            merged += 1
        print(("elapse: ", time.time() - _s))
    print(("add %d entries in the %s['%s']." % (count, database, collection)))
    print(("merge %d entries in the %s['%s']." % (merged, database, collection)))
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
                                 ("lig6", pymongo.ASCENDING),
                                 ("geo_flag", pymongo.ASCENDING),
                                 ("ss_flag", pymongo.ASCENDING)
                                 ])
    if not merged == 0:
        print("=====WARNING====")
        print(("Duplicate complexes(%d) occure in the active learning mode. Should never happen." % merged))
    dump_databse(database_name=database,
                 outpath=outpath,
                 user=user, pwd=pwd)


def push_models(model, model_dict, database, collection,
                user=False, pwd=False,
                host="localhost", port=27017,
                auth=False,
                outpath='/userarchive/db_backup'):
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    this_model = modelMongo(model=model, **model_dict)
    model_id = {"predictor": this_model.predictor, "len_tot": this_model.len_tot, "constraints": this_model.constraints}
    if not query_one(db, collection, constraints=model_id) == None:
        print("A model of has already existed.")
        print(("force_push?", this_model.force_push))
        if this_model.force_push:
            print("force_push is truned on. pushing...")
            db[collection].insert_one(this_model.document)
    else:
        print("pushing...")
        db[collection].insert_one(this_model.document)
    db[collection].create_index([("predictor", pymongo.ASCENDING),
                                 ("len_tot", pymongo.ASCENDING)
                                 ])
    dump_databse(database_name=database,
                 outpath=outpath,
                 user=user, pwd=pwd)


def push_models_actlearn(step, model, database, collection,
                         user=False, pwd=False,
                         host="localhost", port=27017,
                         auth=False,
                         outpath='/userarchive/db_backup'):
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    actlearn_model = modelActLearn(step=step, model=model)
    if not query_one(db, collection, constraints={"step": step}) == None:
        print("A model of step %d has already existed.")
    else:
        db[collection].insert_one(actlearn_model.document)
    dump_databse(database_name=database,
                 outpath=outpath,
                 user=user, pwd=pwd)


def push_models_published(model, model_dict, database, collection,
                          user=False, pwd=False,
                          host="localhost", port=27017,
                          auth=False,
                          outpath='/userarchive/db_backup'):
    db = connect2db(user, pwd, host, port, database, auth)
    ensure_collection(db, collection)
    this_model = modelPublished(model=model, **model_dict)
    model_id = {"target": this_model.target, "doi": this_model.doi, }
    if not query_one(db, collection, constraints=model_id) == None:
        print("A model of has already existed. Please double check...")
    else:
        print("pushing...")
        db[collection].insert_one(this_model.document)
    db[collection].create_index([("dio", pymongo.ASCENDING),
                                 ("target", pymongo.ASCENDING)
                                 ])
    dump_databse(database_name=database,
                 outpath=outpath,
                 user=user, pwd=pwd)
