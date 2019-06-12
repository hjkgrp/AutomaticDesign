import pickle
import numpy as np
from datetime import datetime

mongo_attr_from_run_undef = ["name", "metal", "ox", "spin", "eqlig", "axlig1", "axlig2",
                             "lig1", "lig2", "lig3", "lig4", "lig5", "lig6",
                             "alpha", "functional", "basis", "status", 'converged', 'charge',
                             'terachem_version']
mongo_attr_from_run_nan = ["energy", "ss_target", "ss_actual", 'alphaHOMO', 'betaHOMO',
                           'alphaLUMO', 'betaLUMO']
mongo_attr_other = ["date", "geo_type", "geo_flag", "ss_flag" "opt_geo", "init_geo", "prog_geo", "RACs", "ligcharge"]
mongo_attr_id = ["metal", "ox", "spin", "eqlig", "axlig1", "axlig2",
                 "lig1", "lig2", "lig3", "lig4", "lig5", "lig6",
                 "alpha", "functional", "basis", 'converged']  ### keys that identify a complex in matching.


class tmcMongo():
    def __init__(self, tag, this_run=False, document=False, geo_type=False):
        if not this_run:
            if document:
                try:
                    this_run = pickle.loads(document["dftrun"])
                except:
                    raise ValueError("The input document does not contain a DFTrun object.")
            else:
                raise ValueError("Either a DFTrun object or a tmcMongo object is required as an input.")
        self.document = {}
        self.id_doc = {}
        self.tag = tag
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
        this_run.get_descriptor_vector()
        descriptor_dict = {}
        try:
            for ii, ele in enumerate(this_run.descriptor_names):
                descriptor_dict.update(({ele: this_run.descriptors[ii]}))
        except:
            pass
        self.RACs = descriptor_dict
        self.construct_document()
        self.construct_identity()

    def deserialize_dftrun(self):
        return pickle.loads(self.dftrun)

    def construct_document(self):
        for attr in mongo_attr_from_run_undef + mongo_attr_from_run_nan + mongo_attr_other:
            self.document.update({attr: getattr(self, attr)})
        self.document.update({"dftrun": self.dftrun})
        self.document.update({"tag": self.tag})
        self.document.update({"date": self.date})
        self.document.update({"RACs": self.RACs})

    def construct_identity(self):
        for attr in mongo_attr_id:
            self.id_doc.update({attr: getattr(self, attr)})

    def back_to_mAD(self):
        pass
