import os
import pickle
import numpy as np

attr_actlearn = ["model", "step"]


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


class modelActLearn(MLModel):

    def __init__(self, step, model=False, document=False):
        MLModel.__init__(self, model=model, document=document)
        self.step = step

    def construct_document(self):
        for attr in attr_actlearn:
            try:
                self.document.update({attr: getattr(self, attr)})
            except:
                pass
