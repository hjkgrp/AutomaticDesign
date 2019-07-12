import os
import json
from molSimplifyAD.utils.pymongo_tools import push2db


def deserialize_inputs(filein):
    with open(filein, "r") as fo:
        args_dict = json.load(fo)
    print("db push setup: ", args_dict)
    return args_dict


def check_dbinputs(args_dict):
    keys_required = ["database", "collection", "tag", "subtag", "auth"]
    publication = False
    if set(keys_required) <= set(args_dict.keys()):
        for key in keys_required:
            globals().update({key: args_dict[key]})
    else:
        raise KeyError("missing necessary keys to push data. Required inputs are: ", keys_required)
    if "publication" in args_dict.keys():
        publication = args_dict["publication"]
    return database, collection, tag, subtag, auth, publication


def push_run(args):
    args_dict = deserialize_inputs(args.push)
    database, collection, tag, subtag, auth, publication = check_dbinputs(args_dict)
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
            auth=auth, all_runs_pickle=False)
