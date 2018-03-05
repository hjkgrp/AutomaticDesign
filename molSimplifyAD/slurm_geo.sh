#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres gpu:1
#SBATCH --time 40:00:00
#SBATCH -A p-che140105
#SBATCH -p normal



export OMP_NUM_THREADS=1
module load terachem
export TeraChem="/cstor/xsede/users/xs-jpjanet"

fullpath="$1"


echo "full path is $fullpath"
generalpath=`echo $(dirname $fullpath) | sed "s,/*[^/]\*$,,"`
#echo "gen path is $generalpath"

gennumpath=$(basename $generalpath)
#echo "gen path is $generalpath"

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
opt_geo_path=$generalpath/optimized_geo/$gennumpath/$namebase.xyz
prog_geo_path=$generalpath/prog_geo/$gennumpath/$namebase.xyz
initial_geo_path=$generalpath/initial_geo/$gennumpath/$namebase.xyz
outpath=$generalpath/geo_outfiles/$gennumpath/$namebase.out
completepath=$generalpath/completejobs/$gennumpath/$namebase.done
scrpath=$generalpath/scr/sp/$gennumpath/
echo "scr will be copied to  $scrpath"
echo "paths set"

echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
mkdir -p scr/geo/$gennumpath
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "outpath is $outpath"
echo "optGpath is $opt_geo_path"

#echo "this current env: $SGE_JOB_SPOOL_DIR"
wf_guess_flag=0
##begin geo-optimization

coordfile=$initial_geo_path
echo $coordfile
if [ -e $prog_geo_path ]; then
	echo "restarting from previously optimized geo"
	coordfile=$prog_geo_path #write for continutation in alpha
	guess_opt="$generalpathn/scr/geo/$gennumpath/$namebase/ca0 $generalpathn/scr/geo/$gennumpath/$namebase/cb0"
	wf_guess_flag=1
	echo "Since there is no hand-on guess, and optgeo exists \n"
	echo "I will try to use the scr value. \n"
fi

if [ $wf_guess_flag -eq 0 ]; then ## see if we load in a guess file
	guess_opt="generate"
	echo "wf from scratch"
fi
if [ -e $inpath ]; then
	rm $inpath
fi
echo "copying from $sourcepath to $inpath"
cp $sourcepath $inpath 
sed -i '/min_coordinates cartesian/d' $inpath 
cat >> $inpath <<-EOF
	scf diis+a 
	coordinates $coordfile
	guess $guess_opt
	end
EOF
echo "inpath is $inpath"
echo "Launching geo calc: $namebase"
terachem $inpath >  $localoutpath
msg=$?
if [ -e $scrpath ]; then
	echo "we found scr: $scrpath"
fi
stringtotest="$scrpath/optim.xyz"
cp $localoutpath $outpath
cp -r scr/geo/$namebase $scrpath
echo "Complete"

