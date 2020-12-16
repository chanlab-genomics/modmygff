#!/bin/bash
#
#PBS -N Stri_CCMP2592
#-#PBS -l nodes=1
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:25:00

#PBS -j oe
#PBS -o /gpfs1/homes/s4430291/chanlab-genomics/modmygff/batch/batch_out/Stri_CCMP2592_out.txt

#CHANGE THIS TO YOUR UQ-FACULTY-SCHOOL group name. 
#USE the groups command to find out your exact group name. 
#PBS -A NCMAS-d85

#PBS -l select=1
DATE=$(date +"%d/%m/%Y %H:%M")
echo "time started  "$DATE
echo ------------------------------------------------------
echo -n 'Job is running on the following nodes '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
export TIMEFORMAT="%E sec"

cd $PBS_O_WORKDIR
pwd

module load python
time python3.6 -m modmygff --anno_path /QRISdata/Q2015/ena_genome_submission/Stri_CCMP2592/S.tridacnidorum_CCMP2592_uniprot_annotated.tsv --gff_path /QRISdata/Q2015/ena_genome_submission/Stri_CCMP2592/S.tridacnidorum_CCMP2592.gff --output_path /QRISdata/Q2015/ena_genome_submission/Stri_CCMP2592/S.tridacnidorum_CCMP2592_ext.gff

echo "time finished "$DATE

#Refer to the table of job environment variables section of the RCC PBS Pro User Guide