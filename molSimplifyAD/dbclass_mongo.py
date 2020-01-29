import os
import sys
import getpass
import pickle
import numpy as np
from datetime import datetime

mongo_attr_from_run_undef = ["name", "metal", "ox", "spin", "liglist",
                             "alpha", "functional", "basis", "status", 'converged', 'charge',
                             'terachem_version', "ligcharge", "dynamic_feature", "geo_check_dict",
                             "scrpath", "outpath", "geo_opt", "geo_check_metrics", "geo_check_metrics_prog",
                             'wavefunction', 'molden', 'iscsd']
mongo_attr_flags = ["geo_flag", "ss_flag", "metal_spin_flag"]
mongo_attr_id = ["metal", "ox", "spin", "ligstr", "alpha", "functional", "basis", 'converged',
                 "energy", "geotype"]  ### keys that identify a complex in matching.
csd_attr_id = ["refcode", "metal", "charge", "spin", "functional", "basis", 'converged',
               "energy", "geotype"]  ### keys that identify a complex in matching for csd complexes.
mongo_attr_from_run_nan = ["energy", "ss_target", "ss_act", 'alphaHOMO', 'betaHOMO',
                           'alphaLUMO', 'betaLUMO', 'gap', 'tot_step', 'tot_time',
                           'e_hist', 'e_delta_hist', 'grad_rms_hist',
                           'grad_max_hist', 'displace_rms_hist',
                           'displace_rms_hist', 'displace_rms_hist',
                           'trust_radius_hist', 'step_qual_hist', 'expected_delE_hist',
                           'water_cont', 'thermo_cont',
                           'solvent', 'vertIP', 'vertEA', 'functionalsSP']
SP_keys = ['solvent', 'vertIP', 'vertEA', 'functionalsSP']
mongo_attr_other = ["date", "author", "geotype", "opt_geo", "init_geo", "prog_geo",
                    "RACs", "initRACs", "dftrun", "tag", "subtag", "unique_name",
                    "publication", "ligstr"]
mongo_attr_actlearn = ["step", "is_training", "status_flag", "target", "descriptors", "opt_geo", "init_geo",
                       "ligcharge", "unique_name", "name"]
mongo_not_web = ["dftrun"]
wfn_basepath = '/home/data/wfn/'
dftrun_basepath = '/home/data/dftrun/'
molden_basepath = '/home/data/molden/'


class TMC():
    '''
    Classes that converts between DFTrun and documents in MongoDB.

    Inputs:
        this_run: DFTrun object.
        document: document from MongoDB.
        geotype: type of the geometry of TM complex.
    Note: To successfully initiate a tmcMongo object, either this_run or document is required as an input.

    Key attributes:
        id_doc: a dictionary that tells the unique identity of a TM complex.
        dftrun: a centeralized path to the DFTrun object.
        this_run: DFT object.
    '''

    def __init__(self, this_run=False, document=False, geotype=False):
        if (this_run and document):
            raise ValueError(
                "Confusion. Either a DFTrun object or a tmcMongo object is required as an input. Not both.")
        if not this_run:
            if document:
                try:
                    self.recover_dftrun(document=document)
                except:
                    print(("failing on: ", document["dftrun"]))
                    raise ValueError("The input document cannot recover a DFTrun object.")
            else:
                raise ValueError("Either a DFTrun object or a tmcMongo object is required as an input.")
        else:
            self.this_run = this_run
        if not geotype:
            self.geotype = "oct" if self.this_run.octahedral else "sqpyr"
        else:
            self.geotype = geotype
        self.document = {}
        self.id_doc = {}
        self.update_fields = list()
        self.get_id_keys()
        self.inherit_from_run()
        self.cal_initRAC()
        self.cal_RAC()
        self.construct_identity()
        self.make_unique_name()

    def inherit_from_run(self):
        ## Get attr from dftrun
        for attr in mongo_attr_from_run_undef:
            setattr(self, attr, "undef")
            try:
                setattr(self, attr, getattr(self.this_run, attr))
            except:
                pass
        self.liglist = self.this_run.liglist
        if "liglist_compact" in self.this_run.__dict__:
            self.ligstr = "_".join(self.this_run.liglist_compact)
        else:
            self.ligstr = "_".join(self.this_run.liglist)
        for attr in mongo_attr_from_run_nan:
            setattr(self, attr, np.nan)
            if attr in self.this_run.__dict__:
                try:
                    setattr(self, attr, float(getattr(self.this_run, attr)))
                except:
                    setattr(self, attr, getattr(self.this_run, attr))
        ## Get job flags
        for attr in mongo_attr_flags:
            setattr(self, attr, np.nan)
        for attr in mongo_attr_flags:
            try:
                setattr(self, attr, getattr(self.this_run, attr))
            except:
                pass
        ## Get geometries
        self.opt_geo = self.this_run.mol.returnxyz() if self.this_run.mol else "undef"
        self.init_geo = self.this_run.init_mol.returnxyz() if self.this_run.init_mol else "undef"
        self.prog_geo = self.this_run.progmol.returnxyz() if self.this_run.progmol else "undef"
        ## Get whether is a geometry optimization
        if "geo_opt" in self.this_run.__dict__:
            self.geo_opt = self.this_run.geo_opt
        if "dict_geo_check" in self.this_run.__dict__:
            self.geo_check_metrics = self.this_run.dict_geo_check
        if "dict_geo_check_prog" in self.this_run.__dict__:
            self.geo_check_metrics_prog = self.this_run.dict_geo_check_prog

    def cal_initRAC(self):
        try:
            self.this_run.get_descriptor_vector(useinitgeo=True)
            descriptor_dict = {}
            try:
                for ii, ele in enumerate(self.this_run.descriptor_names):
                    descriptor_dict.update(({ele: self.this_run.descriptors[ii]}))
            except:
                pass
            self.initRACs = descriptor_dict
        except:
            self.initRACs = {}

    def cal_RAC(self):
        try:
            if ("mol" in self.this_run.__dict__ and not self.this_run.status in [1, 2, 8,
                                                                                 9]) or self.this_run.alpha != 20.0:
                self.this_run.get_descriptor_vector(useinitgeo=False)
            elif "init_mol" in self.this_run.__dict__:
                self.this_run.get_descriptor_vector(useinitgeo=True)
            descriptor_dict = {}
            try:
                for ii, ele in enumerate(self.this_run.descriptor_names):
                    descriptor_dict.update(({ele: self.this_run.descriptors[ii]}))
            except:
                pass
            self.RACs = descriptor_dict
        except:
            self.RACs = {}

    def get_id_keys(self):
        if 'iscsd' in self.this_run.__dict__ and self.this_run.iscsd == True:
            self.id_keys = csd_attr_id
            print(("csd complex: ", self.this_run.iscsd, self.this_run.refcode))
            self.refcode = self.this_run.refcode
        else:
            self.id_keys = mongo_attr_id

    def construct_identity(self):
        for attr in self.id_keys:
            self.id_doc.update({attr: getattr(self, attr)})

    def make_unique_name(self):
        name_ele = []
        for key in self.id_keys:
            name_ele.append(key)
            name_ele.append(str(self.id_doc[key]))
        self.unique_name = '_'.join(name_ele)

    def get_update_fields(self, update_fields):
        if update_fields:
            for field in update_fields:
                if field in list(self.document.keys()):
                    self.update_fields.append(field)

    def recover_dftrun(self, document=False):
        if not document:
            dftrun_file = self.document["dftrun"]
        else:
            if "dftrun" in list(document.keys()):
                dftrun_file = document["dftrun"]
            elif "unique_name" in list(document.keys()) and os.path.isdir(dftrun_basepath + document["unique_name"]):
                dftrun_file = dftrun_basepath + document["unique_name"] + '/dftrun.pkl'
            else:
                raise ValueError("Cannot find any DFTrun object binded to this document.")
        if os.path.isfile(dftrun_file):
            if sys.version_info[0] < 3:
                self.this_run = pickle.load(open(dftrun_file, "rb"))
            else:
                self.this_run = pickle.load(open(dftrun_file, "rb"), encoding="latin1")
            for k in document:
                setattr(self.this_run, k, document[k])  # use document values since those are more reliable.
        else:
            raise ValueError("Cannot recover the DFTrun object.")


class tmcMongo(TMC):
    '''
    For the database of DFTRun class.

    Inputs:
        tag: tag of your data. Recommend to use the project name (may related to the name of your paper).
        subtag: Recommend as the name of your mAD folder (considering we may have many mAD folders for each project).
        publication: The paper to which this complex belongs. Format: LastName_Journal_year (e. g., Duan_JCTC_2019)
        this_run: DFTrun object.
        document: document from MongoDB.
        geotype: type of the geometry of TM complex.
        update_fields: fields that will be updated when merging two documents. Default as an empty list.
    Note: To successfully initiate a tmcMongo object, either this_run or document is required as an input.

    Key attributes:
        document: the document to be inserted in MongoDB (contains DFTrub object).
        id_doc: a dictionary that tells the unique identity of a TM complex.
        dftrun: a pickle-dumped DFTrun object.
    '''

    def __init__(self, tag, subtag, publication=False,
                 this_run=False, document=False,
                 geotype=False, update_fields=False):
        TMC.__init__(self, this_run=this_run, document=document, geotype=geotype)
        self.author = getpass.getuser()
        self.tag = tag
        self.subtag = subtag
        self.publication = publication
        self.date = datetime.now()
        self.write_wfn(force=False)
        self.write_molden(force=False)
        self.write_dftrun(force=False)
        self.construct_document()
        self.get_update_fields(update_fields)

    def construct_document(self):
        for attr in mongo_attr_from_run_undef + mongo_attr_from_run_nan + mongo_attr_other + mongo_attr_flags + self.id_keys:
            if attr in self.__dict__:
                self.document.update({attr: getattr(self, attr)})
        for ii, lig in enumerate(self.liglist):
            self.document.update({"lig%d" % (ii + 1): lig})

    def construct_webdoc(self):
        for key in self.document:
            if not key in mongo_not_web:
                try:
                    self.web_doc.update({key: self.document[key]})
                except:
                    pass

    def write_wfn(self, wfn_basepath=wfn_basepath, force=False):
        wfn_path = wfn_basepath + self.unique_name + '/'
        noexist = False
        if (not os.path.isdir(wfn_path)):
            os.makedirs(wfn_path)
            noexist = True
        if self.this_run.wavefunction:
            for key in self.this_run.wavefunction:
                if self.this_run.wavefunction[key]:
                    if noexist or force:
                        with open(wfn_path + key, "wb") as fo:
                            fo.write(self.this_run.wavefunction[key])
                    self.this_run.wavefunction.update({key: wfn_path + key})
                    self.wavefunction.update({key: wfn_path + key})

    def write_molden(self, molden_basepath=molden_basepath, force=False):
        molden_path = molden_basepath + self.unique_name + '/'
        noexist = False
        if 'molden' in self.this_run.__dict__ and not self.this_run.molden == 'undef':
            if (not os.path.isdir(molden_path)):
                os.makedirs(molden_path)
                noexist = True
            if noexist or force:
                if self.this_run.molden:
                    with open(molden_path + "molden.molden", "w") as fo:
                        fo.write(self.this_run.molden)
            self.this_run.molden = molden_path + "molden.molden"
            self.molden = molden_path + "molden.molden"

    def write_dftrun(self, dftrun_basepath=dftrun_basepath, force=False):
        dftrun_path = dftrun_basepath + self.unique_name + '/'
        noexist = False
        if (not os.path.isdir(dftrun_path)):
            os.makedirs(dftrun_path)
            noexist = True
        if noexist or force:
            with open(dftrun_path + "dftrun.pkl", "wb") as fo:
                try:
                    pickle.dump(self.this_run, fo, protocol=2)
                except KeyboardInterrupt:
                    pickle.dump(self.this_run, fo, protocol=2)
                    fo.flush()
        self.dftrun = dftrun_path + "dftrun.pkl"


class tmcActLearn(TMC):
    '''
    For active learning database.

    Inputs:
        step: the step of active learning process.
        is_training: used in train, validation or test
        status_flag: status of the record. 0 as good to use.
        target: target for the learning.
        descriptors: features used in ML predicting target. Default as the RACs for the input DFTrun object.
        this_run: DFTrun object.
        document: document from MongoDB.
        geotype: type of the geometry of TM complex.
        update_fields: fields that will be updated when merging two documents. Default as an empty list.
    Note: To successfully initiate a tmcMongo object, either this_run or document is required as an input.

    Key attributes:
        document: the document to be inserted in MongoDB (contains DFTrub object).
        id_doc: a dictionary that tells the unique identity of a TM complex.
    '''

    def __init__(self, step, is_training, status_flag,
                 target=np.nan, descriptors=False,
                 this_run=False, document=False,
                 geotype=False, update_fields=False):
        TMC.__init__(self, this_run=this_run, document=document, geotype=geotype)
        self.step = step
        self.is_training = is_training
        self.status_flag = status_flag
        self.target = target
        if not descriptors:
            self.descriptors = self.RACs
        else:
            self.descriptors = descriptors
        self.construct_document()
        self.get_update_fields(update_fields)

    def construct_document(self):
        for attr in self.id_keys + mongo_attr_actlearn + mongo_attr_flags:
            if attr in self.__dict__:
                self.document.update({attr: getattr(self, attr)})
        for ii, lig in enumerate(self.liglist):
            self.document.update({"lig%d" % (ii + 1): lig})
