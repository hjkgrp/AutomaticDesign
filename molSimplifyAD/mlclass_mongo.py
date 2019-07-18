import os
import pickle
import numpy as np
from datetime import datetime
import getpass

attr_actlearn = ["model", "step", ]
attr_mongo_in = ["predictor", "hyperparams", "score_train", "score_test",
                 "target_train", "target_test", "pred_train", "pred_test",
                 "len_train", "len_test", "len_tot"]
attr_mongo_undef = ["model", "date", "author"]


class MLModel():

    def __init__(self, model=False, document=False):
        if (model and document):
            raise ValueError(
                "Confusion. Either a model object or a db document is required as an input. Not both.")
        if not model:
            if document:
                try:
                    self.recover_model(document=document)
                except:
                    raise ValueError("The input document cannot recover a Keras model.")
            else:
                raise ValueError("Either a model object or a db document is required as an input.")
        else:
            self.loaded_model = model
        self.model = pickle.dumps(self.loaded_model)
        self.document = {}

    def recover_model(self, document):
        self.loaded_model = pickle.loads(document["model"])


class modelMongo(MLModel):

    def __init__(self, model=False, document=False, **kwargs):
        MLModel.__init__(self, model=model, document=document)
        for attr in attr_mongo_in:
            try:
                setattr(self, attr, kwargs[attr])
            except:
                print("Error. Key %s does not exist in the input dictionary." % attr)
        if "force_push" in kwargs.keys() and kwargs["force_push"]:
            self.force_push = True
        else:
            self.force_push = False
        self.author = getpass.getuser()
        self.date = datetime.now()
        self.construct_document()

    def construct_document(self):
        for attr in attr_mongo_in + attr_mongo_undef:
            try:
                self.document.update({attr: getattr(self, attr)})
            except:
                pass


class modelActLearn(MLModel):

    def __init__(self, step, model=False, document=False):
        MLModel.__init__(self, model=model, document=document)
        self.step = step
        self.construct_document()

    def construct_document(self):
        for attr in attr_actlearn:
            try:
                self.document.update({attr: getattr(self, attr)})
            except:
                pass
