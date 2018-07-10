
import freesasa
import openbabel

class DerivedClassifierT(freesasa.Classifier):
    def classify(self, residueName, atomName):
        print(atomName)
        if atomName.strip() == "CR" or  atomName.strip() == "FE" or  atomName.strip() == "MN":
            print('metal')
            return 'metal'
        elif atomName.strip() == "CO":
            return 'metal'
        else:
            return 'Not-metal'
    def radius(self, residueName, atomName):
        if atomName.strip() == "N": # Nitrogenre.match('\s*N',atomName)
            return 1.55
        if atomName.strip() == "C": # Carbon
            return 1.70
        if atomName.strip() == "O": # Oxygen
            return 1.52
        if atomName.strip() == "S": # Sulfur
            return 1.80
        if atomName.strip() == "P": # Phosphorous
            return 1.80
        if atomName.strip() == "CL": # CHLORINE
            return 1.75
        if atomName.strip() == "H": # Hydrogen
            return 1.10
        if atomName.strip() == "FE": # Iron
            return 2.05
        if atomName.strip() == "CR": # Cr
            print('cr found')
            return 2.05
        if atomName.strip() == "MN": # Mn
            return 2.05
        if atomName.strip() == "CO": # Co
            return 2.00
        if atomName.strip() == "NI": # Ni
            return 2.00
        return 0.0

def get_area(this_run,basename):
    
    path_dictionary=setup_paths()
    outpath = path_dictionary["pdb_path"] + basename + '.pdb'
    
    # convert to pdb
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("xyz")
    obConversion.SetOutFormat("pdb")
    OBMol = openbabel.OBMol()
    obConversion.ReadFile(OBMol, this_run.init_geopath)
    obConversion.WriteFile(OBMol, outpath)

    # measure free SA
    dc = DerivedClassifierT()
    myopt = {'halt-at-unknown': False,
     'hetatm': True,
     'hydrogen': True,
     'join-models': False,
     'skip-unknown': False}
    structure = freesasa.Structure(outpath,classifier = dc, options = myopt)
    structure.setRadiiWithClassifier(dc)

    result = freesasa.calc(structure).totalArea()
    this_run.area = result

