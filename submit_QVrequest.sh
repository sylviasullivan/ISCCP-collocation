#!/bin/sh
#!/rigel/home/scs2229/.conda/envs/mcscfl/bin/python
##
#SBATCH --account=glab                  # The account name for the job.
#SBATCH --job-name=QVdaily              # The job name
#SBATCH -c 12                           # The number of cpu cores to use.
#SBATCH --time=3:00:00                  # The time the job will take to run.
#SBATCH --mem=30gb
#SBATCH --mail-type=NONE                # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=scs2229@columbia.edu

source activate ncplot
python /rigel/home/scs2229/top-secret/MCS_clim/ERAint_QVRetrieve.generated.py


# End of script

