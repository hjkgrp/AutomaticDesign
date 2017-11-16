#$ -S bin/bash		#shell
#$ -l h_rt=24:00:00	#runtime max
#$ -l h_rss=8G		#memory req
#$ -q gpus	        #gpus
#$ -l gpus=1            #
#$ -pe smp 1 		#number parrallel jobs
#  -fin *.py
#  -fin *.sh
#  -fout scr    	# copy back to current dir
#  -fout outfiles
##################CHOOSE A METAL ######################
fullpath="$1"

echo "full path is $fullpath"
generalpath=`echo $(dirname $fullpath) | sed "s,/*[^/]\*$,,"`
echo "gen path is $generalpath"

gennumpath=$(basename $generalpath)
echo "gen path is $generalpath"

generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`
echo "gen path is $generalpath"

generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`
echo "gen path is $generalpath"

namebase=`echo $fullpath | sed "s/[.]in//"| sed "s:.*/::"`

echo "Begining Gibraltar geometry optimization run"
echo "gen path is $generalpath"
export OMP_NUM_THREADS=1
module load cuda
module load terachem
echo "gen path  $generalpath "
echo " gen num path $gennumpath"
echo "base name is  $namebase "
sourcepath=$fullpath
inpath=$generalpath/infiles/$gennumpath/$namebase.in
opt_geo_path=$generalpath/optimized_geo/$gennumpath/$namebase.xyz
initial_geo_path=$generalpath/initial_geo/$gennumpath/$namebase.xyz
outpath=$generalpath/outfiles/$gennumpath/$namebase.out
completepath=$generalpath/completejobs/$gennumpath/$namebase.done
scrpath=$generalpath
echo "scr will be copied to  $scrpath"
echo "paths set"
## sett local, workdir related paths and copy files
localoutpath=$SGE_O_WORKDIR/$namebase.out
localinpath=$SGE_O_WORKDIR/$namebase.in
coordfile=$SGE_O_WORKDIR/$namebase.xyz
cp $initial_geo_path $coordfile 
echo "copied"
echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
mkdir -p scr/sp
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "this current env: $SGE_JOB_SPOOL_DIR"
echo "this SGE WORKDIR: $SGE_O_WORKDIR"

wf_guess_flag=0
##begin geo-optimization

if [ $wf_guess_flag -eq 0 ]; then ## see if we load in a guess file
	guess_opt="generate"
	echo "wf from scratch"
fi
cp $sourcepath $inpath 
cat >> $inpath <<-EOF
	coordinates $coordfile
	guess $guess_opt
	end
EOF
## copy infile back to path
cp $inpath $localinpath

echo "Launching SP calc: $namebase"

terachem $localinpath >  $localoutpath
echo "Complete"

## copy back complete cases 
cp $localinpath /home/jp/tfac/issue.issue
cp $localoutpath $outpath
cp -r scr $scrpath
