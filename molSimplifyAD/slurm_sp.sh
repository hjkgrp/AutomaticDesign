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
echo "gen sp path is $generalpath"

namebase=`echo $fullpath | sed "s/[.]in//"| sed "s:.*/::"`

echo "Begining sp calcualtion run"
echo "general sp path is $generalpath"

echo "gen path = $generalpath" 
echo "namebase = $namebase" 
echo "TeraChem basis dir = $TeraChem"
sourcepath=$fullpath
inpath=$generalpath/sp_infiles/$gennumpath/$namebase.in
outpath=$generalpath/sp_outfiles/$gennumpath/$namebase.out
scrpath=$generalpath/scr/sp/$gennumpath/
localoutpath=$namebase.out
echo "sp scr will be copied to  $scrpath"
echo "sp paths set"

echo "sp inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
mkdir -p scr/sp/$gennumpath/
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "outpath is $outpath"
echo "optGpath is $opt_geo_path"
echo "daemon is in"
pwd

echo "Launching sp calc: $namebase"
$TeraChem/bin/terachem $inpath >  $outpath
mv $localoutpath $outpath
mv scr/sp/$gennumpath/$namebase $scrpath
echo "Complete sp"

