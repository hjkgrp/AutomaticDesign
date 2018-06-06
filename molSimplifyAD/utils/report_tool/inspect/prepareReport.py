import numpy as np
from string import template
import os, sys, shutil, glob, subprocess
from molSimplify.Classes.mol3D import *
from molSimplify.Classes.ligand import *
from molSimplify.Scripts.geometry import *
from pymol.cgo import *
from pymol import cmd
from pymol import preset
from pymol import util

def loadMols(initialPath, finalPath, reportPath):
    # load in geos
    finalMol = mol3D()
    initialMol = mol3D()

    initialMol.readfromxyz(initialPath + '/initial.xyz')
    finalMol.readfromxyz(finalPath + '/final.xyz')

    if not os.path.isdir(reportPath):
        os.makedirs(reportPath)
    filesToCopy = ['main','header']
    for files in filesToCopy:
        shutil.copy(texSource + files + '.tex', reportPath + files + '.tex')
    
    return initialMol, finalMol





def processInitialAndFinalToXyz(initialPath, initialMol, finalPath, finalMol):
    # center on metal
    initialMol.translate([-i for i in initialMol.getAtom(initialMol.findMetal()[0]).coords()])
    finalMol.translate([-i for i in finalMol.getAtom(finalMol.findMetal()[0]).coords()])

    print('######## INITIAL ########')
    initialMolOrient, theta, u = orient(initialMol)
    print('######## FINAL ########')
    finalMolOrient = rotate_around_axis(finalMol,[0,0,0],u[0],-1*theta[0])
    finalMolOrient = rotate_around_axis(finalMol,[0,0,0],u[1],-1*theta[1])
    finalMolOrient = rotate_around_axis(finalMol,[0,0,0],u[2],-1*theta[2])
    print('--reorientation done--')

    initialMolOrient.writexyz(initialPath + '/initialOrient.xyz')
    finalMolOrient.writexyz(finalPath + '/finalOrient.xyz')


# In[44]:


def orient(mol):
    new_mol, theta1, u1 = rot(mol, 2)
    new_mol, theta2, u2 = rot(new_mol, 1)
    new_mol, theta3, u3 = rot(new_mol, 0)
    return new_mol, [theta1, theta2, theta3], [u1, u2, u3]


# In[45]:


def rot(mol, component):
    # component is 0,1,2 for x,y,z
    print('------Component: ' + str(component) + '------')
    maxComponent = 0
    Vect = []
    cInd = 0

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
    print('\nafter rot the conatom is located at ' + str(VectNew))       
    print(VectNew) 
    print('\n')

    return new_mol, theta, u
    


# In[54]:


def draw3d(load_path, save_path):
    cmd.delete('all')

    ## generate axes
    w = 0.06 # cylinder width 
    l = 3 # cylinder length
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
    util.cba
    preset.ball_and_stick(selection='all', mode=1)
    cmd.zoom('all', 1000, 0, 0)
    cmd.set('depth_cue', 0)

    # https://www.andre-gaschler.com/rotationconverter/
    view = [0.8535534, -0.3535534,  0.3826834,       0.5205981,  0.5497281, -0.6532815,       0.0205981,  0.7568349,  0.6532815,         0.000000000,    0.000000000,  -22.395387650,         0.000147820,    0.000066757,    0.000249624,        18.378971100,   26.411804199,  -25.000000000 ]

    cmd.set_view(view)
    cmd.png(save_path + '1.png', 500, 500, 300, 0, 0)
    
    #view2 from the other side
    for i in [2,3,4,5,6,7,8]:
        view[i] = view[i]*(-1)
    
    cmd.set_view(view)
    cmd.png(save_path + '2.png', 500, 500, 300, 0, 0)



# In[78]:


def compileTex():
    print('\n######## COMPILE ########')
    os.chdir('reportName')
    print('Compile in folder:')
    print(os.getcwd())
    cmd_str  = 'pdflatex main.tex'
    p_sub = subprocess.Popen(cmd_str, shell = True, stdout = subprocess.PIPE)
    os.chdir('..')


# In[79]:


def generateReport(initialPath, finalPath, reportPath):
    # paths
    initialPath = initialPath
    finalPath = finalPath
    reportPath = reportPath + "/"

    initialMol, finalMol = loadMols(initialPath, finalPath, reportPath)
    processInitialAndFinalToXyz(initialPath, initialMol, finalPath, finalMol)

    draw3d(initialPath + '/initialOrient.xyz', initialPath + '/initial_render')
    draw3d(finalPath + '/finalOrient.xyz', finalPath + '/final_render')

    compileTex()


# In[80]:


generateReport("initial", "final", "reportName")

