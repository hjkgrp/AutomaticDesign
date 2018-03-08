#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres gpu:1
#SBATCH --time 40:00:00
#SBATCH -A p-che140073
#SBATCH -p normal



export OMP_NUM_THREADS=1
module load terachem
export TeraChem="/cstor/xsede/users/xs-jpjanet"

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
outpath=$generalpath/sp_outfiles/$gennumpath/$namebase.out
scrpath=$generalpath/scr/sp/$gennumpath/
echo "scr will be copied to  $scrpath"
echo "paths set"

echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
mkdir -p scr/sp/
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "outpath is $outpath"
echo "optGpath is $opt_geo_path"
echo "daemon is in"
pwd

localoutpath=$namebase.out

echo "Launching geo calc: $namebase"
terachem $inpath >  $localoutpath
msg=$?
mv $localoutpath $outpath
mv scr/sp/$namebase $scrpath
echo "Complete"

