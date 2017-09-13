#$ -S bin/bash		#shell
#$ -cwd			#return job from cwd
#$ -l h_rt=24:00:00	#runtime max
#$ -l h_rss=8G		#memory req
#$ -q gpus	        #gpus
#$ -l gpus=1            #
#$ -pe smp 1 		#number parrallel jobs
#  -fin *.py
#  -fin *.sh
#  -fout scr    	# copy back to current dir
#  -fout infiles
#  -fout outfiles
#  -fout optimized_geo
#  -fout completejobs
##################CHOOSE A METAL ######################
fullpath="$1"

echo $fullpath
generalpath=`echo $(dirname $fullpath) | sed "s,/*[^/]\*$,,"`
gennumpath=$(basename $generalpath)
generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`
generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`

namebase=`echo $fullpath | sed "s/[.]in//"| sed "s:.*/::"`

echo "Begining Gibraltar geometry optimization run"
export OMP_NUM_THREADS=1
module load cuda
module load terachem
#outpath="$SGE_O_WORKDIR/outfiles/"
echo $generalpath 
echo $gennumpath
echo $namebase 
sourcepath=$fullpath
inpath=$generalpath/infiles/$gennumpath/$namebase.in
opt_geo_path=$generalpath/optimized_geo/$gennumpath/$namebase.xyz
initial_geo_path=$generalpath/initial_geo/$gennumpath/$namebase.xyz
outpath=$generalpath/outfiles/$gennumpath/$namebase.out
completepath=$generalpath/completejobs/$gennumpath/$namebase.done

echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
mkdir -p scr/sp
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
#echo "this current env: $SGE_JOB_SPOOL_DIR"
wf_guess_flag=0
##begin geo-optimization

coordfile=$initial_geo_path
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
echo "Launching SP calc: $namebase"
terachem $inpath >  $outpath
echo "Complete"

