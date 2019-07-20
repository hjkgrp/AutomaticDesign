import os
import json
from molSimplifyAD.utils.pymongo_tools import push2db, push_complex_actlearn


def deserialize_inputs(filein):
    with open(filein, "r") as fo:
        args_dict = json.load(fo)
    print("db push setup: ", args_dict)
    return args_dict


def check_dbinputs(args_dict,get_step=None):
    db_keyword_dict = {}
    keys_required = ["database", "collection", "tag", "subtag", "auth"]
    if get_step != None:
        keys_required += ["step"]
        step = args_dict["step"]
    else:
        step = False
    publication = False
    update_fields = False
    dict_vars = keys_required
    if set(keys_required) <= set(args_dict.keys()):
        for key in keys_required:
            locals().update({key: args_dict[key]})
    else:
        raise KeyError("missing necessary keys to push data. Required inputs are: ", keys_required)
    dict_vars += ["publication","update_fields"]
    if "publication" in args_dict.keys() and (args_dict["publication"] != False):
        publication = args_dict["publication"]
    if "update_fields" in args_dict.keys() and (args_dict["update_fields"] != False):
        update_fields = args_dict["update_fields"]
    for var in dict_vars:
        db_keyword_dict.update({str(var): locals()[var]})
    return db_keyword_dict
    #return database, collection, tag, subtag, auth, publication, update_fields


def push_run(args):
    args_dict = deserialize_inputs(args.push)
    db_keyword_dict = check_dbinputs(args_dict)
    database = db_keyword_dict['database']
    collection = db_keyword_dict['collection']
    tag = db_keyword_dict['tag']
    subtag = db_keyword_dict['subtag']
    auth = db_keyword_dict['collection']
    publication = db_keyword_dict['publication']
    update_fields = db_keyword_dict['update_fields']
    #print(database)
    #database, collection, tag, subtag, auth, publication, update_fields = check_dbinputs(args_dict)
    if auth:
        if (args.user and args.pwd):
            user = args.user
            pwd = args.pwd
        else:
            raise KeyError("missing necessary keys to push data. Required inputs are: [user, pwd].")
    else:
        user = False
        pwd = False
    push2db(database=database, collection=collection,
            tag=tag, subtag=subtag, publication=publication,
            user=user, pwd=pwd, host="localhost", port=27017,
            auth=auth, all_runs_pickle=False,
            update_fields=update_fields)

def push_run_act_learn(args):
    args_dict = deserialize_inputs(args.push_act_learn)
    db_keyword_dict = check_dbinputs(args_dict,True)
    database = db_keyword_dict['database']
    collection = db_keyword_dict['collection']
    tag = db_keyword_dict['tag']
    subtag = db_keyword_dict['subtag']
    auth = db_keyword_dict['collection']
    publication = db_keyword_dict['publication']
    update_fields = db_keyword_dict['update_fields']
    step = db_keyword_dict['step']
    #database, collection, tag, subtag, auth, publication, update_fields = check_dbinputs(args_dict)
    if auth:
        if (args.user and args.pwd):
            user = args.user
            pwd = args.pwd
        else:
            raise KeyError("missing necessary keys to push data. Required inputs are: [user, pwd].")
    else:
        user = False
        pwd = False
    from molSimplifyAD.ga_check_jobs import check_all_current_convergence
    _, all_runs, act_learn_dict_lists = check_all_current_convergence(post_all=True)
    ### act_learn_dict is list with 2 elements. First one is oxo, second one is hat
    ordering = ['oxo','hat'] #overwriting collection
    push_complex_actlearn(step=step, all_complexes=act_learn_dict_lists[0], database=database, collection="act_learn_oxo",
                          user=user, pwd=pwd, host="localhost", port=27017,
                          auth=auth)
    push_complex_actlearn(step=step, all_complexes=act_learn_dict_lists[1], database=database, collection="act_learn_hat",
                          user=user, pwd=pwd, host="localhost", port=27017,
                          auth=auth)
