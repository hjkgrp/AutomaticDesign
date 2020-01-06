from molSimplify.Classes.mol3D import mol3D, atom3D
from molSimplify.Scripts.geometry import *
from molSimplifyAD.ga_tools import *
import numpy as np
import os, sys, shutil, time

if len(sys.argv) > 1:
    dry_run = False
    print('warning, dry run is OFF. 5 second sleep engaged (not too late to cancel...) ')
    time.sleep(5)
else:
    dry_run = True
    print('dry run is ON. Files are safe.  2 second sleep engaged...' )
    time.sleep(2) 

joblist = '/home/nandy/intermediate/jobs/hydroxyl_to_add.txt'
submitted_job_dictionary = find_submitted_jobs()
## live jobs:
live_job_dictionary = find_live_jobs()
## conv'd jobs
converged_jobs = find_converged_job_dictionary()
metal_spin_dictionary = spin_dictionary()
metals_list = get_metals()
for job in list(converged_jobs.keys()):
    path_dictionary = setup_paths()
    gene, gen, _, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, _, ahf, basename,_= translate_job_name(job)
    if int(ahf)!= 20:
        continue
    metal = metals_list[metal]
    flag = True
    if spin not in metal_spin_dictionary[metal][ox]:
        print(('skipped '+str(basename)+' because spin '+str(spin)+' not in '+str(metal)+'('+str(ox)+') dictionary.\n'))
        flag = False
        continue
    if int(ahf) == 20 and axlig2 == 'oxo' and flag == True: #Only generate the hydroxo if the oxo converged
        #At this point, advance the path to text parse...
        path_dictionary =advance_paths(path_dictionary,gen)
        hydroxyl_ox = int(ox) - 1
        for spin_value in metal_spin_dictionary[metal][ox]:
            gene_list = basename.split('_')
            gene_list[-1] = str(spin_value)
            current_gene = "_".join(gene_list)
            if not os.path.exists(path_dictionary["geo_out_path"]+current_gene+'.out'):
                print((current_gene+'.out does not exist. Continuing'))
                continue
            with open(path_dictionary["geo_out_path"]+current_gene+'.out') as f:
                data = f.readlines()
                found_conv = False
                for i, lines in enumerate(data):
                    if (str(lines).find('Optimization Converged.') != -1) or (str(lines).find('Converged!') != -1):
                        found_conv = True
                geo = basename+'.xyz'
                modgeo1 = basename.split('_') #Going to make a geo for +/- 1 spin away from oxo val (checking if it is valid)
                modgeo2 = basename.split('_')
                ligs = get_ligands()
                for i, item in enumerate(ligs):
                    if 'hydroxyl' in item:
                        value = str(i)
                modgeo1[5], modgeo1[8] = str(hydroxyl_ox), str(value)
                modgeo2[5], modgeo2[8] = str(hydroxyl_ox), str(value)
                #If the oxo has spin val of 5, spin val of 6 and 4 hydroxyls will be generated.
                upperspin = int(modgeo1[-1])+1
                lowerspin = int(modgeo1[-1])-1
                modgeo1[-1] = str(upperspin) #get the +1 spin case
                modgeo2[-1] = str(lowerspin) #get the -1 spin case 
                modgeo1 = '_'.join(modgeo1)+'.xyz'
                modgeo2 = '_'.join(modgeo2)+'.xyz'
                if found_conv and os.path.exists(path_dictionary["optimial_geo_path"]+geo):
                    print(('Geo is '+str(geo)))
                    print(('Mod geos are '+str(modgeo1)+' '+str(modgeo2)))
                    mymol = mol3D()
                    mymol.readfromxyz(path_dictionary["optimial_geo_path"]+geo)
                    metalval = mymol.findMetal()
                    bondedatoms = mymol.getBondedAtomsSmart(metalval)
                    oxo = mymol.getAtom(-1)
                    oxo_coord = oxo.coords()
                    metal_coord = mymol.getAtom(metalval[0]).coords()
                    bond1_coord = mymol.getAtom(bondedatoms[0]).coords()
                    bond2_coord = mymol.getAtom(bondedatoms[1]).coords()
                    oxo, dxyz = setPdistance(oxo, oxo_coord, metal_coord,1.84)
                    metaloxo = np.array(oxo_coord)-np.array(metal_coord)
                    extra = getPointu(oxo_coord, 1.2, metaloxo)
                    moveup = np.array(extra)-metal_coord
                    val = midpt(bond1_coord,bond2_coord)
                    movevect = np.array(normalize(np.array(val)+moveup-np.array(oxo_coord)))+oxo_coord
                    p = list(movevect)
                    newH = atom3D('H', p)
                    mymol.addAtom(newH)
                    hydrogen = mymol.getAtom(-1)
                    hydrogen, dxyz = setPdistance(hydrogen, hydrogen.coords(), oxo.coords(),1)
                    if not dry_run:
                        if not os.path.exists(path_dictionary["initial_geo_path"]+modgeo1) and (upperspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                            mymol.writexyz(path_dictionary["initial_geo_path"]+modgeo1)
                        elif (upperspin not in metal_spin_dictionary[metal][hydroxyl_ox]):
                            print((str(upperspin)+' spin not in dictionary for '+str(metal)+'('+str(ox)+')'))
                        else:
                            print((str(modgeo1)+' already exists!'))
                        if not os.path.exists(path_dictionary["initial_geo_path"]+modgeo2) and (lowerspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                            mymol.writexyz(path_dictionary["initial_geo_path"]+modgeo2)
                        elif (lowerspin not in metal_spin_dictionary[metal][hydroxyl_ox]):
                            print((str(upperspin)+' spin not in dictionary for '+str(metal)+'('+str(ox)+')'))
                        else:
                            print((str(modgeo2)+' already exists!'))
                    else:
                        if not os.path.exists(path_dictionary["initial_geo_path"]+modgeo1) and (upperspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                            print(('Would have written geom '+str(modgeo1)))
                        elif (lowerspin not in metal_spin_dictionary[metal][hydroxyl_ox]):
                            print((str(lowerspin)+' spin not in dictionary for '+str(metal)+'('+str(ox)+')'))
                        else:
                            print((str(modgeo1)+' already exists!'))
                        if not os.path.exists(path_dictionary["initial_geo_path"]+modgeo2) and (lowerspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                            print(('Would have written geom '+str(modgeo2)))
                        elif (lowerspin not in metal_spin_dictionary[metal][hydroxyl_ox]):
                            print((str(lowerspin)+' spin not in dictionary for '+str(metal)+'('+str(ox)+')'))
                        else:
                            print((str(modgeo2)+' already exists!'))
                infile1 = modgeo1.strip('.xyz')+'.in'
                infile2 = modgeo2.strip('.xyz')+'.in'
                print(('UPPER INFILE EXISTS:',os.path.exists(path_dictionary["job_path"]+infile1)))
                print(('LOWER INFILE EXISTS:',os.path.exists(path_dictionary["job_path"]+infile2)))
                print(('UPPER IN METAL SPIN DICT:',upperspin in metal_spin_dictionary[metal][hydroxyl_ox]))
                print(('LOWER IN METAL SPIN DICT:',lowerspin in metal_spin_dictionary[metal][hydroxyl_ox]))
                #if not (os.path.exists(path_dictionary["job_path"]+infile1) and os.path.exists(path_dictionary["initial_geo_path"]+modgeo1)) and (upperspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                if not (os.path.exists(path_dictionary["job_path"]+infile1)) and (upperspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                    if not os.path.exists(joblist):
                        if not dry_run:
                            with open(joblist,'w') as listf:
                                listf.write(str(path_dictionary["job_path"]+infile1)+'\n')
                        else:
                            print(('Would have added '+str(path_dictionary["job_path"]+infile1)))
                    else:
                        if not dry_run:
                            with open(joblist,'a') as listf:
                                listf.write(str(path_dictionary["job_path"]+infile1)+'\n')
                        else:
                            print(('Would have added '+str(path_dictionary["job_path"]+infile1)))
                    if not dry_run:
                        with open(path_dictionary["job_path"]+basename+'.in','r') as sourcef:
                            sourcelines = sourcef.readlines()
                            with open(path_dictionary["job_path"]+infile1,'w') as newf:
                                for line in sourcelines:
                                    if not ("scr" in line) and (not "end" in line) and (not "spinmult" in line):
                                        newf.write(line)
                                    elif "spinmult" in line:
                                        newf.write("spinmult "+infile1.strip('.in').split('_')[-1]+'\n')
                                    elif "method" in line:
                                        if int(upperspin) == 1:
                                            newf.write("method b3lyp\n")
                                        else:
                                            newf.write(line)
                                newf.write('scrdir scr/geo/gen_0/'+infile1.strip('.in')+'\n')
                        newf.close()
                        sourcef.close()
                    else:
                        print('------------------WOULD HAVE WRITTEN-------------------')
                        with open(path_dictionary["job_path"]+basename+'.in','r') as sourcef:
                            sourcelines = sourcef.readlines()
                            if True:
                            #with open(path_dictionary["job_path"]+infile1,'w') as newf:    
                                for line in sourcelines:
                                    if not ("scr" in line) and (not "end" in line) and (not "spinmult" in line):
                                        print(line)
                                    elif "spinmult" in line:
                                        print(("spinmult "+infile1.strip('.in').split('_')[-1]+'\n'))
                                    elif "method" in line:                                                                             
                                        if int(upperspin) == 1:                                              
                                            print("method b3lyp\n")                                                               
                                        else:                                                                                          
                                            print(line)                                                                           
                                print(('scrdir scr/geo/gen_0/'+infile1.strip('.in')+'\n'))                                           
                        sourcef.close()
                        print('-------------------------------------------------------------')
                #if not (os.path.exists(path_dictionary["job_path"]+infile2) and os.path.exists(path_dictionary["initial_geo_path"]+modgeo2)) and (lowerspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                if not (os.path.exists(path_dictionary["job_path"]+infile2)) and (lowerspin in metal_spin_dictionary[metal][hydroxyl_ox]):
                    if not os.path.exists(joblist):
                        if not dry_run:
                            with open(joblist,'w') as listf:
                                listf.write(str(path_dictionary["job_path"]+infile2)+'\n')
                        else:
                            print(('Would have added '+str(path_dictionary["job_path"]+infile2)))
                    else:
                        if not dry_run:
                            with open(joblist,'a') as listf:
                                listf.write(str(path_dictionary["job_path"]+infile2)+'\n')
                        else:
                            print(('Would have added '+str(path_dictionary["job_path"]+infile2)))
                    if not dry_run:
                        with open(path_dictionary["job_path"]+basename+'.in','r') as sourcef:
                            sourcelines = sourcef.readlines()
                            with open(path_dictionary["job_path"]+infile2,'w') as newf:
                                for line in sourcelines:
                                    if not ("scr" in line) and (not "end" in line) and (not "spinmult" in line):
                                        newf.write(line)
                                    elif "spinmult" in line:
                                        newf.write("spinmult "+infile2.strip('.in').split('_')[-1]+'\n')
                                    elif "method" in line:
                                        if int(lowerspin) == 1:
                                            newf.write("method b3lyp\n")
                                        else:
                                            newf.write(line)
                                newf.write('scrdir scr/geo/gen_0/'+infile2.strip('.in')+'\n')
                        newf.close()
                        sourcef.close()
                    else:
                        print('------------------WOULD HAVE WRITTEN-------------------')
                        with open(path_dictionary["job_path"]+basename+'.in','r') as sourcef:
                            sourcelines = sourcef.readlines()
                            if True:
                            #with open(path_dictionary["job_path"]+infile2,'w') as newf:    
                                for line in sourcelines:
                                    if not ("scr" in line) and (not "end" in line) and (not "spinmult" in line):
                                        print(line)
                                    elif "spinmult" in line:
                                        print(("spinmult "+infile2.strip('.in').split('_')[-1]+'\n'))
                                    elif "method" in line:                                                                             
                                        if int(lowerspin) == 1:                                              
                                            print("method b3lyp\n")                                                               
                                        else:                                                                                          
                                            print(line)                                                                           
                                print(('scrdir scr/geo/gen_0/'+infile2.strip('.in')+'\n'))                                           
                        sourcef.close()
                        print('-------------------------------------------------------------')
