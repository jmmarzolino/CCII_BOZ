
#!/bin/bash -l
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00:00
#SBATCH --array=3-3%1
#SBATCH --output=merge_duplicates
#SBATCH --job-name="Merge"
#SBATCH --partition=batch

module load samtools
cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/BAM
outfile=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f1  | tail -n1`
infile1=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f2  | tail -n1`
infile2=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f3  | tail -n1`
infile3=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f4  | tail -n1`
infile4=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f5  | tail -n1`

samtools merge -@10 $outfile $infile1 $infile2 $infile3 $infile4
samtools index -c $outfile
