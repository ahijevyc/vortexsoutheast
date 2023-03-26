#!/bin/bash

# Run a stormtime (good for getting NARR wind shear calcs done in parallel)
# usage:
# run_cmds.bash Fay 2002090600

N=$1$2
echo $1 $2 > $1$2.txt

cat <<EOS |qsub
#!/bin/csh
#PBS -N $N
#PBS -A NMMM0021
#PBS -q casper@casper-pbs
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=1:mem=1GB,walltime=00:15:00

conda activate 
python /glade/scratch/ahijevyc/vortexsoutheast/scripts/NARR_composite.py $1$2.txt
EOS
