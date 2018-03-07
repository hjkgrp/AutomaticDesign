#$ -S bin/bash		#shell
#$ -l h_rt=24:00:00	#runtime max
#$ -l h_rss=8G		#memory req
#$ -q gpus	        #gpus
#$ -l gpus=1            #
#$ -pe smp 1 		#number parrallel jobs

module load cuda
module load terachem


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
export OMP_NUM_THREADS=1
echo "gen path  $generalpath "
echo " gen num path $gennumpath"
echo "base name is  $namebase "
sourcepath=$fullpath
inpath=$generalpath/infiles/$gennumpath/$namebase.in
outpath=$generalpath/outfiles/$gennumpath/$namebase.out
scrpath=$generalpath/scr/sp/$gennumpath/
echo "scr will be copied to  $scrpath"
echo "paths set"
## set local, workdir related paths and copy files
localoutpath=$namebase.out
localinpath=$namebase.in
mkdir -p scr
mkdir -p scr/sp
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "this current SGE_JOB_SPOOL_DIR: $SGE_JOB_SPOOL_DIR"
echo "this SGE WORKDIR: $SGE_O_WORKDIR"
echo "this SGE_O_PATH: $SGE_O_PATH"
echo "this SGE shell current DIR: $PWD"
wf_guess_flag=0

##begin sp calculation
cp $inpath $localinpath
echo "Launching SP calc: $namebase"
terachem $localinpath >  $localoutpath
echo "Complete"
## copy back complete cases 
mv $localoutpath $outpath
mv -r scr/sp/$namebase $scrpath
