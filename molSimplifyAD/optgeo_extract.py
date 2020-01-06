import sys
import os


def parse_xyz(f):
    s = f.read()
    ss = s.splitlines()
    lines = [line.split() for line in ss]
    el = list()
    x = list()
    y = list()
    z = list()
    controlflag = 0
    for i in lines:
        if controlflag:
            controlflag = 0
            el = list()
            x = list()
            y = list()
            z = list()
        else:
            try:
                z.append(i[3])
                el.append(i[0])
                x.append(i[1])
                y.append(i[2])
            except:
                controlflag = 1
    natoms = len(el)
    crds = list()
    elements = list()
    for i in range(natoms):
        crds.append([x[i], y[i], z[i]])
        elements.append(el[i])
    ret_dict = dict()
    ret_dict['coords'] = crds
    ret_dict['elements'] = elements
    return ret_dict


def extract_file_check(filepath, targetpath):
    return_status = 1  ## incomplete
    if os.path.exists(filepath):
        filename = os.path.basename(filepath)
        name, ext = os.path.splitext(filename)
        # print("read from " + filepath + " to " + targetpath)
        with open(filepath, 'r') as f:
            ret_dict = parse_xyz(f)
        ## analyzingif the file is empty:
        number_of_elements = len(ret_dict['elements'])
        if number_of_elements < 1:
            print(("File " + filepath + ' is empty'))
            return_status = 2  ## empty file
        else:
            with open(targetpath, 'w') as f:
                f.write(str(number_of_elements) + '\n')
                f.write('# extracted from Terachem optimization\n')
                for i, elements in enumerate(ret_dict['elements']):
                    writebuffer = [elements] + ret_dict['coords'][i]
                    f.write(' '.join(s for s in writebuffer) + '\n')
            return_status = 0  ## success
    else:
        print(('Error. File ' + filepath + ' does not exist. Aborting import.'))
    return (return_status)
