#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres gpu:1
#SBATCH --time 40:00:00
#SBATCH -A p-che140073
#SBATCH -p normal


module load  GCC/4.9.2-binutils-2.25
module load binutils/2.25
module load intel
module load MPICH2
module load CUDA/8.0.44
module unload intel
module load intel/2016
export TeraChem="/cstor/xsede/projects/p-che140073/opt/terchem/06292018"
export LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH


export OMP_NUM_THREADS=1
fullpath="$1"


echo "full path is $fullpath"
generalpath=`echo $(dirname $fullpath) | sed "s,/*[^/]\*$,,"`
#echo "gen path is $generalpath"

gennumpath=$(basename $generalpath)
#echo "gennumpath is $gennumpath"

generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`
#echo "gen path is $generalpath"

generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`
echo "gen path is $generalpath"

namebase=`echo $fullpath | sed "s/[.]in//"| sed "s:.*/::"`


echo "Begining calcualtion run"
echo "general path is $generalpath"

echo "gen path = $generalpath" 
echo "namebase = $namebase" 
echo "TeraChem basis dir = $TeraChem"
sourcepath=$fullpath
inpath=$generalpath/infiles/$gennumpath/$namebase.in
outpath=$generalpath/geo_outfiles/$gennumpath/$namebase.out
scrpath=$generalpath/scr/geo/$gennumpath/
localoutpath=$namebase.out
echo "scr will be copied to  $scrpath"
echo "paths set"

echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
mkdir -p scr/geo/$gennumpath/
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "outpath is $outpath"
echo "optGpath is $opt_geo_path"
echo "daemon is in "
pwd

#echo "this current env: $SGE_JOB_SPOOL_DIR"
wf_guess_flag=0
##begin geo-optimization

echo "inpath is $inpath"
echo "Launching geo calc: $namebase"
$TeraChem/bin/terachem $inpath >  $outpath
mv $localoutpath $outpath
mv scr/geo/$gennumpath/$namebase $scrpath
echo "Complete"

