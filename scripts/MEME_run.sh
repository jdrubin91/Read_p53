#PBS -q long2gb
#PBS -S /bin/bash
#PBS -N MEME
#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l mem=10gb

### Set your expected walltime
#PBS -l walltime=100:00:00

### Setting to mail when the job is complete
#PBS -e /scratch/Users/joru1876/qsub_errors/
#PBS -o /scratch/Users/joru1876/qsub_stdo/  

### Set your email address
#PBS -m ae
#PBS -M joru1876@colorado.edu


module load meme_4.10.1_4

### now call your program

/opt/meme/4.10.1_4/bin/meme $fastafile -oc $output -maxsize 10000000 -dna

