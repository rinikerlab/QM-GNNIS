source ../setup.sh

model=/cluster/work/igc/kpaul/projects/small_molecule_multisolvent/MachineLearning/trained_models/ProductionRun_seed_1612_49_ckpt.pt

for solvent in water chloroform methanol dmso benzene thf acetonitrile pyridine acetone; do
    sbatch --array=0-256 -n 2 --tmp=2000 --time=24:00:00 --mem-per-cpu=16000 \
    --output=/cluster/project/igc/kpaul/Metoxy_logs/Methoxy_orca_min_hessian_no_mpi_smd_B3LYP_TZVP/slurm-%A_%a.out \
    --wrap="python run_minimisation_id_orca_min_hessian_no_mpi_B3LYP_TZVP.py --smiles 'COCCO' --solvent $solvent --output /cluster/project/igc/kpaul/Methoxy_results/Methoxy_start_orca_min_hessian_no_mpi_SMD_B3LYP_TZVP_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/project/igc/kpaul/COCCO_conformers.sdf"
done

for solvent in diethylether dioxane; do
    sbatch --array=0-256 -n 2 --tmp=2000 --time=24:00:00 --mem-per-cpu=16000 \
    --output=/cluster/project/igc/kpaul/Metoxy_logs/Methoxy_orca_min_hessian_no_mpi_smd_B3LYP_TZVP/slurm-%A_%a.out \
    --wrap="python run_minimisation_id_orca_min_hessian_no_mpi_B3LYP_TZVP.py --smiles 'COCCO' --solvent $solvent --output /cluster/project/igc/kpaul/Methoxy_results/Methoxy_start_orca_min_hessian_no_mpi_SMD_B3LYP_TZVP_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/project/igc/kpaul/COCCO_conformers.sdf"
done

for solvent in water chloroform methanol dmso hexane; do
    sbatch --array=0-256 -n 2 --tmp=2000 --time=24:00:00 --mem-per-cpu=16000 \
    --output=/cluster/project/igc/kpaul/I4_logs/I4_orca_min_hessian_no_mpi_smd_B3LYP_TZVP/slurm-%A_%a.out \
    --wrap="python run_minimisation_id_orca_min_hessian_no_mpi_B3LYP_TZVP.py --smiles 'COCCOC' --solvent $solvent --output /cluster/project/igc/kpaul/I4_results/I4_start_orca_min_hessian_no_mpi_SMD_B3LYP_TZVP_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/project/igc/kpaul/COCCOC_conformers.sdf"
done