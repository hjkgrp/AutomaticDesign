import argparse
import copy
import glob
from molSimplifyAD.dbclass_mongo import tmcMongo
from molSimplifyAD.utils.pymongo_tools import connect2db, insert, count_find
from molSimplifyAD.ga_tools import get_mulliken


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-user')
    parser.add_argument('-pwd')
    args = parser.parse_args()
    return args


def main():
    args = arg_parser()
    if args.user == None or args.pwd == None:
        raise KeyError("Please use the format python update_db_documents.py -user <username> -pwd <password>.")
    constraints = {"publication": "Duan_JCTC_2019", "converged": False}
    update_fields = ["geo_flag", "ss_flag", "metal_spin_flag"]
    database = "tmc"
    collection = "oct"
    user = args.user
    pwd = args.pwd
    db = connect2db(user, pwd,
                    host="localhost", port=27017,
                    database=database, auth=True)
    print("Number of complexes in the collection:", db.oct.count())
    cursor = db[collection].find(constraints)
    tot = count_find(cursor)
    print("Number of complexes to be updated: ", tot)
    cursor = db[collection].find(constraints)
    print("Are you sure to update %s with constraints %s in %s[%s]? (y/n)" % (str(update_fields), str(constraints),
                                                                              database, collection))
    _in = raw_input()
    if not _in == "y":
        print("Quit. Have a nice day.")
        quit()
    count = 0
    confirmed = False
    for _tmcdoc in cursor:
        print("complex: ", _tmcdoc["unique_name"])
        try:
            _this_tmc = tmcMongo(document=_tmcdoc, tag=_tmcdoc["tag"], subtag=_tmcdoc["subtag"],
                                 publication=_tmcdoc["publication"], update_fields=update_fields)

            ####
            ## Simple modification. You have already know what to update and you **don't** want to update dftrun.
            # Change here. e.g. _this_tmc.document["publication"] = xxx
            ####

            ####
            ## Modify both the documents and dftrun.
            _this_run = copy.deepcopy(_this_tmc.this_run)
            print("!!!!", _this_run.metal_spin_flag)
            current_folder = _this_run.scrpath.strip("optim.xyz")
            multiwfnpath = glob.glob(current_folder + "*.molden")
            if len(multiwfnpath) > 0:
                multiwfnpath = multiwfnpath[0]
                mulliken_spin_list = get_mulliken(multiwfnpath, _this_run.spin, _this_run.liglist[-1])
                print(mulliken_spin_list)
                _this_run.net_metal_spin = mulliken_spin_list[0]
                if len(mulliken_spin_list) > 1:
                    _this_run.net_oxygen_spin = mulliken_spin_list[1]
            else:
                print("No molden path found for this run (" + str(jobs) + ")")
            _this_run.get_check_flags()
            new_tmc = tmcMongo(this_run=_this_run, tag=_tmcdoc["tag"], subtag=_tmcdoc["subtag"],
                               publication=_tmcdoc["publication"], update_fields=update_fields)
            ###

            if not confirmed:
                for key in update_fields:
                    print("=======Key======: ", key)
                    print("Current: ", _tmcdoc[key])
                    print("Change to: ", new_tmc.document[key])
                print("Is this expected? (y/n)")
                _in = raw_input()
                if _in == "y":
                    confirmed = True
                else:
                    print("Quit. Have a nice day.")
                    quit()
            __ = insert(db, collection, new_tmc)
            count += 1
            print(" In progress: %d / %d" % (count, tot))
        except:
            print("Something went wrong. Please check on complex: ", _tmcdoc["unique_name"])
    print("You have changed %d documents in %s[%s]" % (tot, database, collection))


if __name__ == '__main__':
    main()
