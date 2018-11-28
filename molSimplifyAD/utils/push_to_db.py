from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import *
import os, random, shutil, pickle, datetime, sys, time 
from molSimplify.Scripts.geometry import *
from molSimplify.Classes.atom3D import *
from molSimplify.Classes.globalvars import globalvars
from molSimplify.Classes.mol3D import*
from molSimplify.Informatics.geo_analyze import *
from molSimplify.Informatics.coulomb_analyze import *
from molSimplify.Informatics.autocorrelation import *
from molSimplifyAD.post_classes import *
from molSimplifyAD.process_scf import *
from molSimplifyAD.ga_tools import *
from molSimplifyAD.ga_check_jobs import *
import jsonpickle
import openbabel
import molSimplifyAD.db_base_generator
import molSimplifyAD.dbclass 

## load the path
engine_path = isKeyword('db_path')
## start the connection
engine = create_engine(engine_path, echo=True)
## get base creator class
molSimplifyAD.db_base_generator.Base().metadata.create_all(engine, checkfirst=True)
## open session
Session = sessionmaker(bind=engine)
session = Session()
## force post all
## TO DO
## get results from final convergence
final_results, all_runs =  check_final_convergence()
## commit to DB
print('Warning: db push is enabled, attempting commit with tag ' + str(isKeyword('tag')))
time.sleep(5)
suc = 0
fail = 0
for run_c in all_runs.values():
    try:
        suc += 1
        molSimplifyAD.dbclass.add_to_db(session,run_c,isKeyword('tag'))
    except:
        fail += 1
        print('db add failed for ' + str(run_c.name))
print('converted ' + str(suc) + ' cases, failed for ' + str(fail))
print('sql-commit in 10 seconds...')
time.sleep(10)
print('using sql credentials ' + engine_path)
## commit
session.commit()
print('commited!')
