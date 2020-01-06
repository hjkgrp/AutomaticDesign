import os, sys, shutil
from operator import itemgetter
ligInds = {}
with open('eightJobs','r') as f:
        dat = f.readlines()
        for lines in dat:
                l = os.path.basename(lines.split(',')[0])
                this_slot = l.split('_')[3]
                l = l.split('_')[6]
                if not l in list(ligInds.keys()):
                        ligInds.update({l:[this_slot]})
                elif l in list(ligInds.keys()):
                        ll=ligInds[l]
                        ll.append(this_slot)
                        ligInds.update({l:ll})
ligSmiles = []
with open('smu-neutral-homo/finalSmiMonodentate.txt','r') as f:
        dat = f.readlines()
        for lines in dat:
                l = lines.split()[0]
                ligSmiles.append(l)


res = []
for i,ind in enumerate(ligInds.keys()):
        ind = int(ind)
        res.append([ind,ligSmiles[ind],len(ligInds[str(ind)]),sorted(ligInds[str(ind)])])
res = sorted(res, key=itemgetter(2),reverse=True)
with open('eightLigandsCharged.txt','w') as f:
        for lines in res:
                f.write(" ".join(str(x) for x in lines) + '\n')
ligInds = {}
with open('goodJobs','r') as f:
        dat = f.readlines()
        for lines in dat:
                l = os.path.basename(lines.split(',')[0])
                this_slot = l.split('_')[3]
                l = l.split('_')[6]
                if not l in list(ligInds.keys()):
                        ligInds.update({l:[this_slot]})
                elif l in list(ligInds.keys()):
                        ll=ligInds[l]
                        ll.append(this_slot)
                        ligInds.update({l:ll})
ligSmiles = []
with open('smu-neutral-homo/finalSmiMonodentate.txt','r') as f:
        dat = f.readlines()
        for lines in dat:
                l = lines.split()[0]
                ligSmiles.append(l)


res = []
for i,ind in enumerate(ligInds.keys()):
        ind = int(ind)
        res.append([ind,ligSmiles[ind],len(ligInds[str(ind)]),sorted(ligInds[str(ind)])])
res = sorted(res, key=itemgetter(2),reverse=True)
with open('goodLigandsCharged.txt','w') as f:
        for lines in res:
                f.write(" ".join(str(x) for x in lines) + '\n')
