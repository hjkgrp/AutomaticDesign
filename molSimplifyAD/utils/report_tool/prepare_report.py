
import numpy as np
from string import Template
import os, sys, shutil, glob, subprocess
#from molSimplifyAD.ga_tools import *
from molSimplifyAD.ga_io_control import *
from molSimplify.Classes.mol3D import *
from molSimplify.Classes.ligand import *
from molSimplify.Scripts.geometry import *
from molSimplify.Informatics.misc_descriptors import *
from molSimplify.Informatics.graph_analyze import *
from molSimplify.Informatics.geo_analyze import *

allowedPymol = False

try:
        from pymol.cgo import *
        from pymol import cmd
        from pymol import preset
        from pymol import util
        allowedPymol = True
except ImportError:
        #raise ImportError('connection to pymol unsuccessful.')
        print('connection to pymol unsuccessful.')

def flProcess(flag_list):
    formatted_list = ""
    print(flag_list)
    
    if flag_list == "None":
        flag_list = ["None"]
    if isinstance(flag_list,str):
        flag_list = [flag_list]
    for fl in flag_list:
        formatted_list  += verbatimize(str(fl))+' '
    return formatted_list

def verbatimize(st):
    return "\\begin{verbatim}" + st + "\\end{verbatim}"
            
def loadMols(initialPath, finalPath):
    # load in geos
    finalMol = mol3D()
    initialMol = mol3D()
    try:
        initialMol.readfromxyz(initialPath)
    except:
        print("can't load initial mol ")
    try:
        finalMol.readfromxyz(finalPath)
    except:
        print("can't load final mol ")

    return initialMol, finalMol


def processInitialAndFinalToXyz(workdir, initialMol, finalMol):
    # center on metal
    initialMol.translate([-i for i in initialMol.getAtom(initialMol.findMetal()[0]).coords()])
    finalMol.translate([-i for i in finalMol.getAtom(finalMol.findMetal()[0]).coords()])

    #print('######## INITIAL ########')
    initialMolOrient, theta, u = orient(initialMol)
    #print('######## FINAL ########')
    finalMolOrient = rotate_around_axis(finalMol,[0,0,0],u[0],-1*theta[0])
    finalMolOrient = rotate_around_axis(finalMol,[0,0,0],u[1],-1*theta[1])
    finalMolOrient = rotate_around_axis(finalMol,[0,0,0],u[2],-1*theta[2])
    #print('--reorientation done--')

    initialMolOrient.writexyz(workdir + 'initial/initialOrient.xyz')
    finalMolOrient.writexyz(workdir + 'final/finalOrient.xyz')


def orient(mol):
    new_mol, theta1, u1 = rot(mol, 2)
    new_mol, theta2, u2 = rot(new_mol, 1)
    new_mol, theta3, u3 = rot(new_mol, 0)
    return new_mol, [theta1, theta2, theta3], [u1, u2, u3]


def rot(mol, component):
    # component is 0,1,2 for x,y,z
    #print('------Component: ' + str(component) + '------')
    maxComponent = 0
    Vect = []
    cInd = 0
    initialLiglist, initialLigdents, initialLigcons = ligand_breakdown(mol)
    for conAtoms in initialLigcons:
        thisComponent = abs(mol.getAtom(conAtoms[0]).coords()[component])
        if maxComponent < thisComponent:
            Vect = mol.getAtom(conAtoms[0]).coords()
            cInd = conAtoms[0]
    
    vectSize = norm(Vect)
    if component == 0:
        target = [vectSize, 0, 0]
    elif component == 1:
        target = [0, vectSize, 0]
    elif component == 2:
        target = [0, 0, vectSize]
    
    thetaold = vecangle(Vect, target)
    
    theta, u = rotation_params(Vect, [0,0,0], target)
    new_mol = rotate_around_axis(mol, [0,0,0], u, (-1) * theta)
    
    VectNew = mol.getAtom(cInd).coords()
    #print('\nafter rot the conatom is located at ' + str(VectNew))       
    #print(VectNew) 
    #print('\n')

    return new_mol, theta, u
    

def draw3d(load_path, save_path):
    cmd.delete('all')

    ## generate axes
    w = 0.06 # cylinder width 
    l = 5 # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter

    # start coords, end coords, w is radius, color1 with grad to color2 in rgb
    obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           CYLINDER, 0.0, 0.0, 0.0,   -l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   -l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   -l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
          ]
    # load axes into picture
    cmd.load_cgo(obj, 'axes')

    # load structure and define properties 
    cmd.load(load_path)
    cmd.bg_color(color="white")
    preset.ball_and_stick(selection='all', mode=1)

    util.cbaw
    cmd.color("gray","symbol c")
    cmd.select('hats', "symbol h")
    cmd.set("sphere_scale",0,'hats')
    cmd.set("stick_h_scale",1)
    cmd.set("stick_radius",0.70)

    cmd.zoom('all', 1000, 0, 0)
    cmd.set('depth_cue', 0)

    # https://www.andre-gaschler.com/rotationconverter/
    view = [0.8535534, -0.3535534,  0.3826834,       0.5205981,  0.5497281, -0.6532815,       0.0205981,  0.7568349,  0.6532815,         0.000000000,    0.000000000,  -22.395387650,         0.000147820,    0.000066757,    0.000249624,        18.378971100,   26.411804199,  -25.000000000 ]

    cmd.set_view(view)
    cmd.zoom('all', 0.5, 0, 1)
    cmd.png(save_path + '1.png', 500, 500, 300, 0, 0)
    
    #view2 from the other side
    for i in [2,3,4,5,6,7,8]:
        view[i] = view[i]*(-1)
    
    cmd.set_view(view)
    cmd.zoom('all', 0.5, 0, 1)
    cmd.png(save_path + '2.png', 500, 500, 300, 0, 0)


def makeReportFiles(workdir, repDict):
        
    # initiate file copy
    filesToCopy = ['main','header']
    for files in filesToCopy:
        shutil.copy(get_texsource_file(files + '.tex'), workdir + files + '.tex')
  
    # template text with repDict
    filein = open(get_texsource_file('summaryTemplate.tex'),'r')
    filein.seek(0)
    src = Template(filein.read())
    
    #do the substitution
    srcd = src.safe_substitute(repDict)
    with open(workdir + 'summary.tex','w') as f:
        f.write(srcd)
    filein.close()

    
    
def compileTex(workdir,reportPath):
    #print('\n######## COMPILE ########')
    cwd = os.getcwd()
    os.chdir(workdir)
    #print('Compile in folder:')
    #print(os.getcwd())
    cmd_str  = 'pdflatex -interaction=nonstopmode main.tex'
    #print(cmd_str)
    f = open("texoutputFile","wb")
    p_sub = subprocess.call(cmd_str,shell=True,stdout=f)
    f.close()
    #print('done, copying to reportPath')
    os.chdir(cwd)
    shutil.copy(workdir + 'main.pdf',reportPath)


def generateReport(initialPath, finalPath, reportPath,customDict=dict(),octahedral=True):
    # get placeholder dictionary
    repDict = basicRepdict()
    # write custom elemenes to repDict
    repDict = dict(repDict, **customDict)
    
    # make paths safely
    workdir = 'tempReport/' # to be replaced
    workdir = makeCleanWorkdir(workdir)

    # load the two geometries
    initialMol, finalMol = loadMols(initialPath, finalPath)

    # get tex info
    repDict = makeGeoReportDictionary(repDict, initialMol, finalMol, octahedral)

    # rotate and align the geos
    try:
        processInitialAndFinalToXyz(workdir, initialMol, finalMol)
    except:
        print('rotation failed')


    
    # make 3D renders
    if allowedPymol:
        try:
            draw3d(workdir + 'initial/initialOrient.xyz', workdir + 'initial/initial_render')
            draw3d(workdir + 'final/finalOrient.xyz', workdir + 'final/final_render')
        except:
            shutil.copy(get_texsource_file('failsafe.png'), workdir +'initial/initial_render1.png')
            shutil.copy(get_texsource_file('failsafe.png'), workdir +'initial/initial_render2.png')
            shutil.copy(get_texsource_file('failsafe.png'), workdir +'final/final_render1.png')
            shutil.copy(get_texsource_file('failsafe.png'), workdir +'final/final_render2.png')
    else:
        print('No pymol found')
        # implement failsafe fails HERE
        shutil.copy(get_texsource_file('failsafe.png'), workdir +'initial/initial_render1.png')
        shutil.copy(get_texsource_file('failsafe.png'), workdir +'initial/initial_render2.png')
        shutil.copy(get_texsource_file('failsafe.png'), workdir +'final/final_render1.png')
        shutil.copy(get_texsource_file('failsafe.png'), workdir +'final/final_render2.png')

    # move the tex files from templates to the report folder
    makeReportFiles(workdir, repDict)
    
    # write the report
    compileTex(workdir,reportPath)
    
    # now remove workdir
    #print('removing '+ workdir)
    shutil.rmtree(workdir)

def makeCleanWorkdir(workdir):
    org_name = workdir
    counter = 0
    while os.path.isdir(workdir):
        #print 'Warning: ' +workdir +' already exists, generating unique key...'
        workdir =  org_name.rstrip('/') +'_'+ str(counter) + '/'
        counter+=1
    ensure_dir(workdir)
    ensure_dir(workdir+'initial/')
    ensure_dir(workdir+'final/')
   
    return workdir

def basicRepdict():
    # this function 
    # makes a blank dictionary 
    # return dictionary
    repDict = dict()
    # update run info from MAD
    repDict.update({"NAME":"someComplex"})
    repDict.update({"METAL":"someMetal"})
    repDict.update({"LIGS":"reallyGreatLigands"})
    repDict.update({"LIGFORMULAES":"reallyGreatLigandFormulae"})
    repDict.update({"OX":"someOxidationState"})
    repDict.update({"SPIN":"someSpin"})
    repDict.update({"STATUS":"someStatus"})
    repDict.update({"HFX":"someHFX"})
    repDict.update({"s2":"somes2"})
    return repDict



    
def makeGeoReportDictionary(repDict,initialMol,finalMol,octahedral):

    # get ligand formula:
    #if True:
    try:
        if octahedral:  # can only work for oct currently
            initial_axnames, initial_eqnames = getLigFormulae(initialMol)
            print(initial_axnames)
            print(initial_eqnames)
            initial_formulae = "/".join([initial_eqnames[0]]+initial_axnames[0:2])
        else:
            initial_formulae = ""
    except:
        initial_formulae = "er"
        
    repDict.update({"LIGFORMULAES":initial_formulae}) 

    
    # measure bonds
    try:

        init_ax_dist, init_eq_dist = getOctBondDistances(initialMol)
        init_av_ax = '{0:.2f}'.format(np.mean(init_ax_dist))
    except:
        init_av_ax = "er"
    try:
        init_av_eq = '{0:.2f}'.format(np.mean(init_eq_dist))
    except:
        init_av_eq = "er"

    repDict.update({"initialBL":"/".join([init_av_eq,init_av_ax])}) 
    
    
    try:

        final_ax_dist, final_eq_dist = getOctBondDistances(finalMol)
        final_av_ax = '{0:.2f}'.format(np.mean(final_ax_dist))
    except:
        final_av_ax = "er"
    try:
        final_av_eq = '{0:.2f}'.format(np.mean(final_eq_dist))
    except:
        final_av_eq = "er"

    repDict.update({"finalBL":"/".join([final_av_eq,final_av_ax])}) 
    
    # get denticity
    try: 
            resdict = generate_all_ligand_misc(initialMol,loud=False)
            init_ax_dent = [resdict["result_ax"][0]]
            init_eq_dent = [resdict["result_eq"][0]]
            print(init_ax_dent)

    except:
            init_ax_dent = ['er']
            init_eq_dent = ['er']

                
    repDict.update({"initialDent":"/".join([ str(i) for i in init_eq_dent+init_ax_dent])})
    try: 
            resdict = generate_all_ligand_misc(finalMol,loud=False)
            final_ax_dent = [resdict["result_ax"][0]]
            final_eq_dent = [resdict["result_eq"][0]]
    except:
            final_ax_dent = ['er']
            final_eq_dent = ['er']


    repDict.update({"finalDent":"/".join([ str(i) for i in final_ax_dent+final_eq_dent])})
    
    # geo check calls
    try:
        if octahedral:
                init_flag_oct, init_flag_list, init_dict_oct_info = initialMol.IsOct()
        else:
                init_flag_oct, init_flag_list, init_dict_oct_info = initialMol.IsStructure()

        init_AD = '{0:.2f}'.format(init_dict_oct_info["oct_angle_devi_max"])

    except:
        init_AD = 'er'
    repDict.update({"initialAD":init_AD})
    try:
    #if True:
        if octahedral:
                final_flag_oct, final_flag_list, final_dict_oct_info = finalMol.IsOct(init_mol=initialMol)
        else:
                final_flag_oct, final_flag_list, final_dict_oct_info = finalMol.IsStructure(init_mol=initialMol)

        final_AD = '{0:.2f}'.format(final_dict_oct_info["oct_angle_devi_max"])

    except:
        final_AD = 'er'
    repDict.update({"finalAD":final_AD})
    
    # fail list 
    try:
    #if True:
        repDict.update({"initialFL":flProcess(init_flag_list)})
    except:
        repDict.update({"initialFL":'er'})
    try:
    #if True:
       repDict.update({"finalFL":flProcess(final_flag_list)})
    except:
       repDict.update({"finalFL":'er'})
    
    

    
    # rmsd 
    
    try:
        rmsd = '{0:.2f}'.format(final_dict_oct_info["rmsd_max"])
        rmsd = rmsd.zfill(1)
    except:
        rmsd = 'er'
    repDict.update({"RMSD":rmsd})
    return repDict


