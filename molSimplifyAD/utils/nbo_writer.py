import pandas as pd
import sys
import os
import shutil
import glob
write = False
csv_file = sys.argv[1]
write_path = sys.argv[2]
if write_path[0] != '/':
    print('need absolute path to write to')
    sard
write_path = write_path.rstrip('/')
print(('FINALWRITEPATH',write_path))
dataset = pd.read_csv(csv_file)
for i, row in dataset.iterrows():
    #print(row['chem_name'])
    #### This gets the chemical name
    adjusted_name = row['chem_name'].replace('ahf_20','s')
    if write:
        ### This makes a folder for each chemname
        if not os.path.exists(write_path+'/'+adjusted_name):
            os.mkdir(write_path+'/'+adjusted_name)
    destination = write_path+'/'+adjusted_name+'/'
    ##### this copies over the optimized geo #####
    if write:
        ### copies the geo over with the name
        shutil.copy(row['geopath'],destination+adjusted_name+'.xyz')
    
    splitpath = row['geopath'].split('/')
    basepath = '/'.join(splitpath[0:4])

    ####### the scr path is a bit harder to identify due to multiple scr files possibly being present #####
    # currently getting geo/gen_0/
    wfnpath = basepath+'/scr/geo/gen_0/'+row['name']+'*'
    possible_wfn_paths = glob.glob(wfnpath)
    if len(possible_wfn_paths) > 1:
        print('***********************')
        for possible_path in possible_wfn_paths:
            if int(row['spin']) == 1:
                if os.path.exists(possible_path+'/c0'):
                    print('singlet path found')
                    final_wfn_path = possible_path
                    break
                else:
                    print('singlet path NOT found')
            else:
                if os.path.exists(possible_path+'/ca0'):
                    print('other spin path found')
                    final_wfn_path = possible_path
                    break
                else:
                    print('other spin path NOT found')
        print((row['chem_name']))
        print((row['spin']))
        print(possible_wfn_paths)
    else:
        final_wfn_path = possible_wfn_paths[0]
    if write:
        ### copies wfn files over
        if int(row['spin']) == 1:
            shutil.copy(final_wfn_path+'/c0',destination+'c0')
            guess_string = destination+'c0'
        else:
            shutil.copy(final_wfn_path+'/ca0',destination+'ca0')
            shutil.copy(final_wfn_path+'/cb0',destination+'cb0')
            guess_string = destination+'ca0'+' '+destination+'cb0' ## need the guess string for infile
    jobpath = basepath+'/jobs/gen_0/'+row['name']+'.in'
    with open(jobpath) as f:
        data = f.readlines()
        if write:
            with open(destination+adjusted_name+'.in','w') as g:
                for i, line in enumerate(data):
                    if 'scr' in line:
                        g.writelines('scrdir ./scr\n')
                    elif 'run' in line:
                        g.writelines('run energy\n')
                    elif 'method' in line:
                        if int(row['spin']) == 1:
                            g.writelines('method b3lyp\n')
                        else:
                            g.writelines('method ub3lyp\n')
                    else:
                        g.writelines(line)
                g.writelines('nbo yes\n')
                g.writelines('guess '+guess_string+'\n')
                g.writelines('coordinates '+destination+adjusted_name+'.xyz\n')
                g.writelines('end\n')
        else:
            print(('not writing infile, but would normally write to '+ destination+adjusted_name+'.in'))

print('done.')
