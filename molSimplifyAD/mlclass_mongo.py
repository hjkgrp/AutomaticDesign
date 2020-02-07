import os
import pickle
import uuid
import numpy as np
from datetime import datetime
import getpass
from keras.models import load_model

attr_actlearn = ["model", "step", ]
attr_mongo_in = ["predictor", "hyperparams", "score_train", "score_test",
                 "target_train", "target_test", "pred_train", "pred_test",
                 "name_train", "name_test",
                 "len_train", "len_test", "len_tot", "constraints",
                 "history", "tag", "features", "target", "hyperopt",
                 "initialize_weight", "hyperopt_step", "load_latest_model",
                 "fix_architecture", "direct_retrain"]
attr_mongo_undef = ["model", "date", "author"]
model_basepath = '/home/data/models/'


class MLModel():

    def __init__(self, model=False, document=False):
        self.recover_from_doc = False
        if (model and document):
            raise ValueError(
                "Confusion. Either a model object or a db document is required as an input. Not both.")
        if not model:
            if document:
                try:
                    self.recover_model(document=document)
                    self.recover_from_doc = True
                    for key in document:
                        setattr(self, str(key), document[key])
                except:
                    raise ValueError("The input document cannot recover a Keras model.")
            else:
                raise ValueError("Either a model object or a db document is required as an input.")
        else:
            self.loaded_model = model
        self.document = {}

    def recover_model(self, document):
        try:
            self.loaded_model = pickle.loads(document["model"])
        except:
            self.loaded_model = load_model(str(document["model"]))


class modelMongo(MLModel):

    def __init__(self, model=False, document=False, **kwargs):
        MLModel.__init__(self, model=model, document=document)
        for attr in attr_mongo_in:
            if not attr in self.__dict__:
                try:
                    setattr(self, attr, kwargs[attr])
                except:
                    if attr in ["predictor", "history", "hyperparams", "constraints", "len_tot", "features", "target"]:
                        raise KeyError("Error: Key %s does not exist in the input dictionary." % attr)
                    else:
                        print(("Warning: Key %s does not exist in the input dictionary." % attr))
        if "force_push" in list(kwargs.keys()) and kwargs["force_push"]:
            self.force_push = True
        else:
            self.force_push = False
        self.author = getpass.getuser()
        self.date = datetime.now()
        self.write_model()
        self.construct_document()

    def construct_document(self):
        for attr in attr_mongo_in + attr_mongo_undef:
            if attr in self.__dict__:
                self.document.update({attr: getattr(self, attr)})
            else:
                print(("warning: %s does not exist." % attr))

    def write_model(self, model_basepath=model_basepath):
        if not self.recover_from_doc:
            created = False
            while not created:
                id = uuid.uuid4()
                model_path = model_basepath + str(id) + '/'
                if not os.path.isdir(model_path):
                    os.makedirs(model_path)
                    created = True
            self.loaded_model.save(model_path + "model.h5")
            self.model = model_path + "model.h5"


class modelActLearn(MLModel):

    def __init__(self, step, model=False, document=False):
        MLModel.__init__(self, model=model, document=document)
        self.step = step
        self.model = pickle.dumps(self.loaded_model)
        self.construct_document()

    def construct_document(self):
        for attr in attr_actlearn:
            if attr in self.__dict__:
                self.document.update({attr: getattr(self, attr)})
            else:
                print(("warning: %s does not exist." % attr))
