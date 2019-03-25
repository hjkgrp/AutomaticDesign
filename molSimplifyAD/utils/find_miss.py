import os, sys, re, glob, subprocess

subprocess.call('grep -H "No CUDA" geo_outfiles/gen_0/*.out > record.txt',shell=True)

list_of_no_cuda = {}
gpu_dict = {}
with open('record.txt','r') as f:
    for l in f.readlines():
        temp_name = l.strip().split(':')[0].split('/')[-1].strip('.out')
        this_match = glob.glob('queue_output/*'+temp_name+'*.o')
        if len(this_match) >0:
            mtime = os.path.getmtime(this_match[0])
            with open(this_match[0],'r') as q:
                for line in q.readlines():
                    if line.find('releasing') != -1:
                        base_line= line.strip()
                        node = base_line.split('-')[1].split()[0]
                        if len(base_line.split(': ')) > 1:
                                gpu = base_line.split(': ')[1]
                        else:
                                gpu='no-gpu'
                    if line.find('Creating') != -1:
                        base_line = line.strip()
                        job_num = base_line.split('/')[-1]
            if not job_num in list_of_no_cuda.keys():
                list_of_no_cuda.update({job_num:[node, gpu,str(mtime),os.getcwd() +'/' + this_match[0]]})
            node_gpu = node+'_'+gpu
            if not node_gpu in gpu_dict.keys():
                gpu_dict.update({node_gpu:1})
            else:
                gpu_dict[node_gpu] += 1
with open('bad_gpu_jobs2.txt','w') as b:
    for keys in list_of_no_cuda.keys():
        tempvar = list_of_no_cuda[keys]
        b.write(str(keys)+','+','.join(tempvar)+'\n')
print(gpu_dict)
