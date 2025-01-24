
source ../setup.sh
model=/cluster/work/igc/kpaul/projects/small_molecule_multisolvent/MachineLearning/trained_models/ProductionRun_seed_1612_49_ckpt.pt

jobidChloroform=$(sbatch --parsable --array=0-53 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent Chloroform --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_Chloroform_cluster_center.h5")

sleep 180

jobidacetone=$(sbatch --parsable --array=0-56 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent acetone --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_acetone_cluster_center.h5")

sleep 180

jobidacetonitrile=$(sbatch --parsable --array=0-60 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent acetonitrile --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_acetonitrile_cluster_center.h5")

sleep 180

jobidEthylacetate=$(sbatch --parsable --array=0-57 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent Ethylacetate --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_Ethylacetate_cluster_center.h5")

sleep 180

jobidTHF=$(sbatch --parsable --array=0-51 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent THF --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_THF_cluster_center.h5")

sleep 180

jobidDCM=$(sbatch --parsable --array=0-50 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent DCM --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_DCM_cluster_center.h5")

sleep 180

jobidEthanol=$(sbatch --parsable --array=0-63 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent Ethanol --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_Ethanol_cluster_center.h5")

sleep 180

jobidMethanol=$(sbatch --parsable --array=0-61 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent Methanol --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_Methanol_cluster_center.h5")

sleep 180

jobidDMSO=$(sbatch --parsable --array=0-59 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B1_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2N)cc1' --solvent DMSO --output /cluster/project/igc/kpaul/MB_results/B1_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B1/B1_DMSO_cluster_center.h5")

sleep 180

jobid=$(sbatch --parsable -n 1 --time=01:00:00 --wrap="bash submission_scripts/submit_new_orca_hessian_no_mpiB2.sh")
