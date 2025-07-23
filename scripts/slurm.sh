#!/bin/bash
#SBATCH --job-name=GoOrDie
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raffaele.iannuzzi@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/group/iorio/Raffaele/SCDRESP_data/data/results/log/elnet_nopoly_%A_%a.out
#SBATCH --error=/group/iorio/Raffaele/SCDRESP_data/data/results/log/elnet_nopoly_%A_%a.err
#SBATCH --array=1-272%50

####################################################
##### set environment
####################################################

source /home/raffaele.iannuzzi/miniconda3/etc/profile.d/conda.sh
conda activate /home/raffaele.iannuzzi/miniconda3/envs/scanpy_backup_v2/

####################################################

DRUG_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" drugs_list.txt)

echo "Processing drug: $DRUG_NAME"

# Run Python script
python train_elastic_net.py "$DRUG_NAME" "training_sets_dict_nopoly.pkl"
                
