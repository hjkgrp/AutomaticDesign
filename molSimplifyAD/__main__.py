import tensorflow as tf
import keras
from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_create_new import *
from molSimplifyAD.ga_resume_run import *
from molSimplifyAD.ga_database_io import push_run, push_run_act_learn
from molSimplifyAD.retrain.retrain_from_db import retrain_and_push
import sys


def main():
    args = sys.argv[1:]
    print(args)
    # print welcome message
    ss = "\n************************************************************"
    ss += "\n***** Welcome to molSimplify-AD ! Let's get started. *****\n"
    ss += "************************************************************\n\n"
    print(ss)

    parser = argparse.ArgumentParser()
    args = parseall(parser)
    print(args)
    args = checkinput(args)

    print(args)

    if args.new:
        print('Starting new run...')
        create_new_run(args)
    elif args.resume:
        print('Resuming existing run...')
        resume_run(args)
    elif args.push:
        print('Push current run to db...')
        push_run(args)
    elif args.push_act_learn:
        print('Pushing active learning info to db...')
        push_run_act_learn(args)
    elif args.retrain:
        predictor = args.retrain
        print(("Model retraining on :", predictor))
        infile = args.infile
        if not args.user or not args.pwd:
            raise KeyError(
                "Please use the format mad -retrain <predictor> -infile <infile> -user <username> -pwd <password> to retrain a model.")
        args_dict = deserialize_json(infile)
        args_dict.update({"predictor": predictor,
                          "user": args.user,
                          "pwd": args.pwd})
        print(("with arguments as: ", args_dict))
        retrain_and_push(args_dict)


if __name__ == '__main__':
    main()
