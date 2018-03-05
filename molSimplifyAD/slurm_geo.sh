#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres gpu:1
#SBATCH --time 40:00:00
#SBATCH -A p-che140105
#SBATCH -p normal
fullpath="$1"


echo $fullpath
namebase=`echo $fullpath | sed "s/[.]in//"| sed "s:.*/::"`
echo "namebase = $namebase" 
generalpathn="/cstor/xsede/users/xs-jpjanet/xstream_redox"
generalpath="/cstor/xsede/users/xs-jpjanet/xstream_redox/"
echo "Beginingi xstream geometry optimization run"
export OMP_NUM_THREADS=1
module load terachem
export TeraChem="/cstor/xsede/users/xs-jpjanet"
echo "gen path = $generalpath" 
echo "namebase = $namebase" 
echo "TeraChem basis dir = $TeraChem"
sourcepath=$generalpathn/jobs/$namebase.in
inpath=$generalpathn/infiles/$namebase.in
opt_geo_path=$generalpathn/optimized_geo/$namebase.xyz
prog_geo_path=$generalpathn/prog_geo/$namebase.xyz

opt_geo_path=`echo "$opt_geo_path" | awk '{print tolower($0)}'`
echo $opt_geo_path

initial_geo_path=$generalpathn/initial_geo/$namebase.xyz
outpath=$generalpathn/geo_outfiles/$namebase.out
localoutpath=$namebase.out
completepath=$generalpathin/completejobs/$namebase.in
scrpath=$generalpathn/scr/geo/$namebase
geoextr=$generalpathn/optgeo_extract.py
echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p $generalpathn/scr
mkdir -p $generalpathn/scr/geo
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
	guess_opt="$generalpathn/scr/geo/$namebase/ca0 $generalpathn/scr/geo/$namebase/cb0"
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
cp -r  scr
echo "Complete"

