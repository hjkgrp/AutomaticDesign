
from sqlalchemy import Column, Integer,Text, Unicode, UnicodeText, String, Numeric, Date, Text
from molSimplifyAD.ga_tools import *
from molSimplifyAD.post_classes import *
from molSimplify.Classes.mol3D import *
from molSimplify.Classes.atom3D import *
from molSimplifyAD.db_base_generator import Base
import jsonpickle
import getpass
from datetime import datetime as dtn
class TMC(Base):
    ## this class ingests a DFT run class
    ## object and populates a holder
    ## class based on the run that can
    ## be written into a sql database
    ## we hold some basic ID parameters
    ## as well as a json of the class
    ## that allows it to be recreated
    ## 
    
    ## name of the database we use:
    __tablename__ = 'TMDATATWO'
    
    ## this is the unique key that will be used by
    ## sql, by default it is row number added (i.e. seq ints)
    id = Column(Integer, primary_key=True)
    
    ## compulsory records per key
    name = Column(Text)
    metal = Column(Text)
    geo = Column(Text)
    ox = Column(Integer)
    spin = Column(Integer)
    axlig1 = Column(Text)
    axlig2 = Column(Text)
    eqlig = Column(Text)
    aHF = Column(Numeric)    
    author = Column(Text)
    tag = Column(Text)
    date = Column(Date)
    
    
    ## optional records
    energy = Column(Text, nullable = True)
    opt_geo = Column(Text,nullable = True)
    init_geo = Column(Text,nullable = True)
    status = Column(Integer,nullable = True)
    geo_flag = Column(Text,nullable = True)
    reconstructor = Column(Text(length=4294967295), nullable = True)
    spin_squ_target  = Column(Numeric, nullable = True)  
    spin_squ_act = Column(Numeric, nullable = True) 
    
    def __init__(self, run_class,tag):
    
        ## convert compulsory labels
        self.name = run_class.name
        self.spin = int(run_class.spin)
        self.ox = int(run_class.ox)
        self.metal = run_class.metal
        self.axlig1 = run_class.axlig1
        self.axlig2 = run_class.axlig2
        self.eqlig = run_class.eqlig
        self.aHF = run_class.alpha
        self.date = dtn.now()
        self.tag = tag
        self.author = getpass.getuser()

        if run_class.octahedral:
            self.geo = "oct"
        else:
            self.geo = "sqpyr"
        
        ## convergence information
        self.status = int(run_class.status)
        self.geo_flag = str(run_class.flag_oct) 
        self.spin_squ_target = float(run_class.ss_target)
        self.spin_squ_act = float(run_class.ss_act)
        
        ## get energy 
        self.energy = float(run_class.energy)
        print('before mapping on save ')
        print(len(str(( jsonpickle.encode(run_class)))))
        ## get json of class for reconstruction
        #run_class.descriptors = []
        #run_class.descriptor_names = []
        #run_class.mol = []
        print('after desc purge on save ')
        print(len(str(( jsonpickle.encode(run_class)))))
        self.reconstructor = jsonpickle.encode(run_class)
        print('\n')
        print('on save ')
        print(len(str((self.reconstructor))))
        print('*****')
        # capture geos
        if run_class.mol: 
            self.opt_geo = run_class.mol.returnxyz()
        if run_class.init_mol: 
            self.init_geo = run_class.init_mol.returnxyz()
    def reinstate_run_class(self):
        print('\n')
        print('length of stored ')
        print(len(str(self.reconstructor)))
        with open('rc-1-loaded.json','w') as f:
            f.write(self.reconstructor)
        print('*****')
        print('length of json ')
        print(len(jsonpickle.decode(self.reconstructor)))
        print('*****')
        print('\n')

        print('\n')
        #run_class = jsonpickle.decode(self.reconstructor)
        run_class = 0

        print('*****')
        return(run_class)
########################
           
def add_to_db(session,run_class,tag):
    ## needs to call session commit after all
    ## changes are done!
    dbc = TMC(run_class,tag=tag)
    session.add(dbc)

########################
    
def query_to_db(session):
    ## query the DB for all TM data
    db_query = session.query(TMC)
    return(db_query)

#######################
    
def match_against_db(session,metal,ox,spin,geo,eqliq,axlig1,axlig2,eqlig,aHF):
    ## check a run potential TMC against the DB
    ## returns a lost of TMC objects
    db_query = session.query(TMC)
    matches = []
    for r in db_query:
        if (r.metal == metal) and    \
           (r.ox == ox) and    \
           (r.spin == spin) and    \
           (r.geo == geo) and    \
           (r.eqlig == eqlig) and    \
           (r.axlig1 == axlig1) and  \
           (r.axlig2 == axlig2) and  \
           (r.aHF == aHF):
            matches.append(r)
    return(matches)