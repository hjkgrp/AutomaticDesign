import argparse
from molSimplifyAD.dbclass_mongo import tmcMongo
from molSimplifyAD.utils.pymongo_tools import connect2db, insert, count_find


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
    constraints = {"publication": "Duan_JCTC_2019"}
    update_fields = ["publication"]
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
            this_tmc = tmcMongo(document=_tmcdoc, tag=_tmcdoc["tag"], subtag=_tmcdoc["subtag"],
                                publication=_tmcdoc["publication"], update_fields=update_fields)

            ####
            # Change here. e.g. _tmc.document["publication"] = xxx
            ####

            if not confirmed:
                for key in update_fields:
                    print("=======Key======: ", key)
                    print("Current: ", _tmcdoc[key])
                    print("Change to: ", this_tmc.document[key])
                    print("Is this expected? (y/n)")
                    _in = raw_input()
                    if _in == "y":
                        confirmed = True
                    else:
                        print("Quit. Have a nice day.")
                        quit()
            __ = insert(db, collection, this_tmc)
            count += 1
            print(" In progress: %d / %d" % (count, tot))
        except:
            print("Something went wrong. Please check on complex: ",  _tmcdoc["unique_name"])
    print("You have changed %d documents in %s[%s]" % (tot, database, collection))


if __name__ == '__main__':
    main()
