import os
import glob
import numpy as np
from molSimplify.job_manager.tools import textfile


def bind_direct(this_run, jobname, basedir, case, keyinout, suffix=''):
    case_attr = case + suffix
    setattr(this_run, case_attr, False)
    outfile = basedir + "/" + "%s_%s" % (jobname, case) + "/" + "%s_%s.out" % (jobname, case)
    if os.path.isfile(outfile):
        setattr(this_run, case_attr, np.nan)
        output = textfile(outfile)
        v = output.wordgrab([keyinout[0]], [keyinout[1]], last_line=True)[0]
        if not v == None:
            setattr(this_run, case_attr, v)


def bind_with_search(this_run, jobname, basedir, case, keyinout):
    setattr(this_run, case, False)
    search_dir = basedir + "/" + "%s_%s" % (jobname, case)
    if os.path.isdir(search_dir):
        setattr(this_run, case, dict())
        for dirpath, dirs, files in os.walk(search_dir):
            for file in sorted(files):
                if file.split('.')[-1] == 'out':
                    outfile = dirpath + '/' + file
                    key = file.strip('out').strip(jobname).strip('.')
                    output = textfile(outfile)
                    energy = output.wordgrab([keyinout[0]], [keyinout[1]], last_line=True)[0]
                    if not energy == None:
                        getattr(this_run, case).update({key: energy})
                    else:
                        getattr(this_run, case).update({key: False})


def bind_water(this_run, jobname, basedir):
    bind_direct(this_run, jobname, basedir,
                case='water',
                keyinout=['C-PCM contribution to final energy:', -2],
                suffix='_cont')


def bind_thermo(this_run, jobname, basedir):
    bind_direct(this_run, jobname, basedir,
                case='thermo',
                keyinout=['Thermal vibrational energy (ZPE + <E>)', -2],
                suffix='_cont')


def bind_solvent(this_run, jobname, basedir):
    bind_with_search(this_run, jobname, basedir,
                     case="solvent",
                     keyinout=['C-PCM contribution to final energy:', -2])


def bind_vertIP(this_run, jobname, basedir):
    this_run.vertIP = False


def bind_ligdissociate(this_run, jobname, basedir):
    bind_with_search(this_run, jobname, basedir,
                     case="dissociation",
                     keyinout=['FINAL ENERGY:', -2])
