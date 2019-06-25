import os
from molSimplifyAD.utils.pymongo_tools import push2db


def deserialize_inputs(filein):
    args_dict = {}
    if os.path.isfile(filein):
        with open(filein, 'r') as fo:
            for line in fo:
                ll = line.split()
                if len(ll) == 2:
                    args_dict.update({str(ll[0]): str(ll[1])})
                else:
                    raise ValueError("length of each line in %s can only contain 2 elements." % filein)
    else:
        raise ValueError("file %s do not exist." % filein)
    print("db push setup: ", args_dict)
    return args_dict


def check_dbinputs(args_dict):
    keys_required = ["database", "collection", "tag", "subtag", "auth"]
    if set(keys_required) <= set(args_dict.keys()):
        for key in keys_required:
            globals().update({key: args_dict[key]})
    else:
        raise KeyError("missing necessary keys to push data. Required inputs are: ", keys_required)
    return database, collection, tag, subtag, auth


def push_run(args):
    args_dict = deserialize_inputs(args.push)
    database, collection, tag, subtag, auth = check_dbinputs(args_dict)
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
            user=user, pwd=pwd, host="localhost", port=27017,
            auth=auth, all_runs_pickle=False,
            tag=tag, subtag=subtag, web="web")
