import os
import getpass
import pickle
import numpy as np
from datetime import datetime

mongo_attr_from_run_undef = ["name", "metal", "ox", "spin", "lig1", "lig2", "lig3", "lig4", "lig5", "lig6",
                             "alpha", "functional", "basis", "status", 'converged', 'charge',
                             'terachem_version']
mongo_attr_from_run_nan = ["energy", "ss_target", "ss_act", 'alphaHOMO', 'betaHOMO',
                           'alphaLUMO', 'betaLUMO']
mongo_attr_other = ["date", "author", "geo_type", "geo_flag", "ss_flag", "opt_geo", "init_geo", "prog_geo",
                    "RACs", "initRACs", "ligcharge", "dftrun", "tag", "subtag", "unique_name"]
mongo_attr_id = ["metal", "ox", "spin", "lig1", "lig2", "lig3", "lig4", "lig5", "lig6",
                 "alpha", "functional", "basis", 'converged']  ### keys that identify a complex in matching.
mongo_not_web = ["dftrun"]


class tmcMongo():

    '''
    Classes that converts between DFTrun and documents in MongoDB.

    Inputs:
        tag: tag of your data. Recommend to use the project name (may related to the name of your paper).
        subtag: Recommend as the name of your mAD folder (considering we may have many mAD folders for each project).
        this_run: DFTrun object.
        document: document from MongoDB.
        geo_type: type of the geometry of TM complex.
    Note: To successfully initiate a tmcMongo object, either this_run or document is required as an input.

    Key attributes:
        document: the document to be inserted in MongoDB (contains DFTrub object).
        web_doc: a simplified document to be inserted in MongoDB for the web interface (DFTrun object excluded).
        id_doc: a dictionary that tells the unique identity of a TM complex.
        dftrun: a pickle-dumped DFTrun object.
    '''

    def __init__(self, tag, subtag, this_run=False, document=False, geo_type=False):
        if (this_run and document):
            raise ValueError(
                "Confusion. Either a DFTrun object or a tmcMongo object is required as an input. Not both.")
        if not this_run:
            if document:
                try:
                    this_run = pickle.loads(document["dftrun"])
                except:
                    raise ValueError("The input document does not contain a DFTrun object.")
            else:
                raise ValueError("Either a DFTrun object or a tmcMongo object is required as an input.")
        self.document = {}
        self.web_doc = {}
        self.id_doc = {}
        self.author = getpass.getuser()
        self.tag = tag
        self.subtag = subtag
        self.date = datetime.now()
        self.dftrun = pickle.dumps(this_run)
        for attr in mongo_attr_from_run_undef:
            setattr(self, attr, "undef")
        for attr in mongo_attr_from_run_undef:
            try:
                setattr(self, attr, getattr(this_run, attr))
            except:
                pass
        for attr in mongo_attr_from_run_nan:
            setattr(self, attr, np.nan)
        for attr in mongo_attr_from_run_nan:
            try:
                setattr(self, attr, float(getattr(this_run, attr)))
            except:
                pass
        if not geo_type:
            self.geo_type = "oct" if this_run.octahedral else "sqpyr"
        else:
            self.geo_type = geo_type
        self.geo_flag = this_run.flag_oct
        try:
            _, _ = float(this_run.ss_act), float(this_run.ss_target)
            self.ss_flag = 1 if (this_run.ss_act - this_run.ss_target) < 1 else 0
            if this_run.spin == 1:
                self.ss_flag = 1
        except TypeError:
            self.ss_flag = "undef"
        try:
            _, _ = int(self.charge), int(self.ox)
            self.ligcharge = int(self.charge) - int(self.ox)
        except ValueError:
            self.ligcharge = "undef"
        self.opt_geo = this_run.mol.returnxyz() if this_run.mol else "undef"
        self.init_geo = this_run.init_mol.returnxyz() if this_run.init_mol else "undef"
        self.prog_geo = this_run.progmol.returnxyz() if this_run.progmol else "undef"
        this_run.get_descriptor_vector(useinitgeo=True)
        descriptor_dict = {}
        try:
            for ii, ele in enumerate(this_run.descriptor_names):
                descriptor_dict.update(({ele: this_run.descriptors[ii]}))
        except:
            pass
        self.initRACs = descriptor_dict
        this_run.get_descriptor_vector(useinitgeo=False)
        descriptor_dict = {}
        try:
            for ii, ele in enumerate(this_run.descriptor_names):
                descriptor_dict.update(({ele: this_run.descriptors[ii]}))
        except:
            pass
        self.RACs = descriptor_dict
        self.construct_identity()
        self.make_unique_name()
        self.construct_document()
        self.construct_webdoc()

    def deserialize_dftrun(self):
        return pickle.loads(self.dftrun)

    def construct_document(self):
        for attr in mongo_attr_from_run_undef + mongo_attr_from_run_nan + mongo_attr_other:
            try:
                self.document.update({attr: getattr(self, attr)})
            except:
                pass

    def construct_webdoc(self):
        for attr in mongo_attr_from_run_undef + mongo_attr_from_run_nan + mongo_attr_other:
            if not attr in mongo_not_web:
                try:
                    self.web_doc.update({attr: getattr(self, attr)})
                except:
                    pass

    def construct_identity(self):
        for attr in mongo_attr_id:
            self.id_doc.update({attr: getattr(self, attr)})

    def make_unique_name(self):
        name_ele = []
        for key in mongo_attr_id:
            name_ele.append(key)
            name_ele.append(str(self.id_doc[key]))
        self.unique_name = '_'.join(name_ele)

    def write_wfn(self, this_run, wfn_basepath='/data/wfn/'):
        wfn_path = wfn_basepath + self.unique_name + '/'
        if not os.path.isdir(wfn_path):
            os.makedirs(wfn_path)
        for key in this_run.wavefunction:
            if this_run.wavefunction[key]:
                with open(wfn_path + key, "wb") as fo:
                    fo.write(this_run.wavefunction[key])
            this_run.wavefunction.update({key: wfn_path + key})
        return this_run

    def back_to_mAD(self):
        pass #TODO
