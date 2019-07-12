import os
import pickle
import numpy as np

attr_actlearn = ["model", "step"]


class MLModel():

    def __init__(self, model):
        self.model = pickle.dumps(model)
        self.document = {}


class modelActLearn(MLModel):

    def __init__(self, model, step):
        MLModel.__init__(self, model=model)
        self.step = step

    def construct_document(self):
        for attr in attr_actlearn:
            try:
                self.document.update({attr: getattr(self, attr)})
            except:
                pass


