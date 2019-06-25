from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_create_new import *
from molSimplifyAD.ga_resume_run import *
from molSimplifyAD.ga_database_io import push_run
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
    
if __name__ == '__main__':
    main()
