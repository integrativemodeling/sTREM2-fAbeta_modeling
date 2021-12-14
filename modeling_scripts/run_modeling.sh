#!/bin/env bash
#$ -cwd             ## use current working directory
#$ -l scratch=20G  ## needs 200 GB of /scratch space
#$ -S /bin/bash
#$ -t 1-9
#$ -l h_rt=336:00:00
#$ -pe smp 8


## 0. In case TMPDIR is not set, e.g. on development nodes, set
##    it to local /scratch, if it exists, otherwise to /tmp
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

parent_dir=`pwd`

## 1. load modules

module load Sali
module load imp/2.13.0-x86_64
module load mpi

## 2. run directory for this independent run

interval=6
i=$(expr $SGE_TASK_ID \* $interval)
RUNDIR='only_xl_'$i
initial=$i
final=$(expr $initial + $interval)

echo $initial
echo $final

## 3. run sampling for this run
if [ ! -d $RUNDIR ]; then
    mkdir $RUNDIR
    cd $RUNDIR    

    current_dir=`pwd`
    
    ## 4. Use a temporary working directory
    cd "$TMPDIR"
    echo "hello"
    
    ## 5. Copy input files from global disk to local scratch
    cp -r $parent_dir/data .
    cp $parent_dir/modeling.py .
    
    ## 6. run the modeling script now
    mpirun -np $NSLOTS python modeling.py $initial $final > only_xl.log

    ## 7. Move output files back to global disk
    mv only_xl* $current_dir
    mv *.log $current_dir
    echo $HOSTNAME
    echo $TMPDIR
    cd $parent_dir

    ## 8. End-of-job summary
    [[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

fi
