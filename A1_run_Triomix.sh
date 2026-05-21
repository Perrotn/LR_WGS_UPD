sample_data="${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_STR_Report/Reference_files/Bioinfo-LR_SampleData.csv"
my_virtualenv="/home/pernic02/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate"

version=0.9

rm -f ${HOME}/scratch/slurm_TrioMix_v${version}.out

sbatch --output=slurm_TrioMix_v${version}.out \
    ${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/A2_slurm_TrioMix.sh \
    --sample_data "$sample_data" \
    --version "$version" \
    --my_virtualenv "$my_virtualenv" \
    --update_only "true"




sample_data="${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_STR_Report/Reference_files/Bioinfo-LR_SampleData.csv"
my_virtualenv="/home/pernic02/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate"

version=0.9.1

rm -f /home/pernic02/scratch/slurm_TrioMix_v${version}.out

sbatch --output=slurm_TrioMix_v${version}.out \
    ${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/A2_slurm_TrioMix.sh \
    --sample_data "$sample_data" \
    --version "$version" \
    --my_virtualenv "$my_virtualenv" \
    --update_only "true"


# grep "HSJ-059-03" /home/pernic02/scratch/slurm_TrioMix_v0.9.1.out
# grep "HSJ-059-03" /home/pernic02/scratch/slurm_TrioMix_v0.9.out
