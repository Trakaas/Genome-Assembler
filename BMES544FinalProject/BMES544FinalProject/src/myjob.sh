#!/bin/bash
#
### tell SGE to use bash for this script
#$ -S /bin/bash
### execute the job from the current working directory, i.e. the directory in which the qsub command is given
#$ -cwd
### join both stdout and stderr into the same file
#$ -j y
### set email address for sending job status
#$ -M cle39@drexel.edu
### project - basically, your research group name with "Grp" replaced by "Prj"
# -P sacanGrp
### request 5 hours of wall clock time "h_rt" = "hard real time" (format is HH:MM:SS, or integer seconds)
#$ -l h_rt=05:00:00
### a hard limit 8 GB of memory per slot - if the job grows beyond this, the job is killed
#$ -l h_vmem=16G
### want nodes with at least 6 GB of free memory per slot
#$ -l m_mem_free=6G
### select the queue all.q, using hostgroup @intelhosts
#$ -q all.q@@intelhosts

. /etc/profile.d/modules.sh

### These four modules must ALWAYS be loaded
module load shared
module load proteus
module load sge/univa
module load gcc

### Whatever modules you used, in addition to the 4 above, 
### when compiling your code (e.g. proteus-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
module load python/3.5-current

export STAGINGDIR=/scratch/$USER/$JOB_ID
mkdir -p $STAGINGDIR

export PROGRAM=/home/cle39/Genome%20Assembler/BMES544FinalProject/BMES544FinalProject/src/BMES544FinalProject.py
cp -R ../* $STAGINGDIR
cd $STAGINGDIR

python3.5 BMES544FinalProject.py -o 10

mv $STAGINGDIR/data/* /home/cle39/Genome%20Assembler/BMES544FinalProject/BMES544FinalProject/data/

/bin/rm -rf $STAGINGDIR
