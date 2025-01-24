
source ../setup.sh
model=/cluster/work/igc/kpaul/projects/small_molecule_multisolvent/MachineLearning/trained_models/ProductionRun_seed_1612_49_ckpt.pt

jobidChloroform=$(sbatch --parsable --array=0-361 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent Chloroform --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_Chloroform_cluster_center.h5")

sleep 180

jobidacetone=$(sbatch --parsable --array=0-447 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent acetone --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_acetone_cluster_center.h5")

sleep 180

jobidacetonitrile=$(sbatch --parsable --array=0-404 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent acetonitrile --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_acetonitrile_cluster_center.h5")

sleep 180

jobidEthylacetate=$(sbatch --parsable --array=0-451 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent Ethylacetate --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_Ethylacetate_cluster_center.h5")

sleep 180

jobidTHF=$(sbatch --parsable --array=0-420 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent THF --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_THF_cluster_center.h5")

sleep 180

jobidDCM=$(sbatch --parsable --array=0-344 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent DCM --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_DCM_cluster_center.h5")

sleep 180

jobidEthanol=$(sbatch --parsable --array=0-427 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent Ethanol --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_Ethanol_cluster_center.h5")

sleep 180

jobidMethanol=$(sbatch --parsable --array=0-431 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent Methanol --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_Methanol_cluster_center.h5")

sleep 180

jobidDMSO=$(sbatch --parsable --array=0-534 -n 2 --tmp=10000 --time=24:00:00 --mem-per-cpu=8000 --output=/cluster/project/igc/kpaul/MB_logs/B3_orca_min_hessian_no_mpi/slurm-%A_%a.out --wrap="python run_minimisation_id_orca_min_hessian_no_mpi.py --smiles 'Fc1ccc(N(C=O)CCc2ccccc2Nc2ccccc(C#N)2)cc1' --solvent DMSO --output /cluster/project/igc/kpaul/MB_results/B3_start_orca_min_hessian_no_mpi_f025 --idx \$SLURM_ARRAY_TASK_ID --model $model --file /cluster/work/igc/kpaul/projects/multisolvent_pub/Simulation/Minimizations/MB_results/B3/B3_DMSO_cluster_center.h5")

sleep 180

jobid=$(sbatch --parsable -n 1 --time=01:00:00 --wrap="bash submission_scripts/submit_new_orca_hessian_no_mpiB4.sh")
