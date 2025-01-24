import argparse

argparser = argparse.ArgumentParser()

argparser.add_argument("--solvent", type=str, required=True, help="Solvent to use")
argparser.add_argument("--model", type=str, required=True, help="Path to the model")
argparser.add_argument("--output", type=str, required=True, help="Output file")
argparser.add_argument("--file", type=str, required=True, help="File with positions")
argparser.add_argument("--idx", type=int, required=True, help="Index of the conformer")
argparser.add_argument(
    "--smiles", type=str, required=True, help="SMILES string of the molecule"
)
args = argparser.parse_args()

### Imports
from MachineLearning.helper_functions import create_gnn_model

from MachineLearning.GNN_Models import GNN3_Multisolvent_embedding_run_multiple_Delta
from implicitml.calculator import (
    ImplicitSolventCalculator,
    Functional,
    BasisSet,
    SolventModel,
    Solvent,
)
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import mdtraj
from rdkit.Geometry import Point3D

import uuid

os.system(f"mkdir -p {args.output}")

# Check if file exists and energies are not infinite
if os.path.isfile(f"{args.output}/{args.solvent}_{args.idx}.npy"):
    print("File exists")
    data = np.load(f"{args.output}/{args.solvent}_{args.idx}.npy", allow_pickle=True)
    if not np.isnan(data[0][1]):
        print("Already calculated")
        exit()
    else:
        print("Recalculating")

solvent_mapping = {
    "water": "tip3p",
    "chloroform": "Chloroform",
    "methanol": "Methanol",
    "dmso": "DMSO",
    "benzene": "Benzol",
    "thf": "THF",
    "acetonitrile": "acetonitrile",
    "pyridine": "pyridine",
    "acetone": "acetone",
    "dcm": "DCM",
    "ethanol": "Ethanol",
    "ethylacetate": "Ethylacetate",
    "diethylether": "Diethylether",
    "dioxane": "14dioxane",
    "hexane": "Hexan",
}

orca_solvent = eval(f"Solvent.{args.solvent.upper()}")
gnn_solvent = solvent_mapping[args.solvent.lower()]

if args.file.endswith(".h5"):
    smiles = args.smiles
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(
        mol,
        numConfs=1,
        randomSeed=42,
        useExpTorsionAnglePrefs=False,
        useBasicKnowledge=True,
    )
    trajectory = mdtraj.load(args.file)
    # Set positions in molecule
    conf = mol.GetConformer(0)
    for i in range(mol.GetNumAtoms()):
        x, y, z = trajectory.xyz[args.idx, i] * 10
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
elif args.file.endswith(".sdf"):
    mol = Chem.SDMolSupplier(args.file, removeHs=False)[args.idx]

results = []


# Create temporary directory manually
tempdir = os.environ.get("TMPDIR", "/tmp")


def create_temp_dir_with_random_name(tempdir):
    uuids = str(uuid.uuid4())
    os.makedirs(f"{tempdir}/{uuids}")
    return f"{tempdir}/{uuids}/"


# add a random uuid to the tempdir
base_tempdir = create_temp_dir_with_random_name(tempdir)

for i in range(mol.GetNumConformers()):
    attempts = 0
    pos = np.zeros((mol.GetNumAtoms(), 3))
    free_energy = np.inf

    smd_orca_pos = np.nan
    smd_free_energy = np.nan
    cosmo_rs_free_energy = np.nan
    gnn_cpcm_free_energy_smd_opt = np.nan
    cpcm_hessian_at_smd_opt = np.nan
    gnn_hessian_at_smd_opt = np.nan
    gnn_smd_free_energy_smd_opt = np.nan
    cpcm_gnn_single_point_at_smd_opt = np.nan
    gnn_opt_pos_0025 = np.nan
    cpcm_gnn_free_energy_gnn_opt = np.nan
    cpcm_hessian_at_gnn_opt = np.nan
    gnn_hessian_at_gnn_opt = np.nan
    cpcm_gnn_single_point_at_gnn_opt = np.nan

    while attempts < 1:
        try:
            model = create_gnn_model(
                mol,
                solvent=gnn_solvent,
                model_class=GNN3_Multisolvent_embedding_run_multiple_Delta,
                model_dict=args.model,
                jit=False,
            )

            new_mol = Chem.Mol(mol)
            new_mol.RemoveAllConformers()
            new_mol.AddConformer(mol.GetConformer(i), assignId=True)
            calc = ImplicitSolventCalculator(new_mol)

            # 1 perform Orca minimization
            print("Optimizing with SMD", flush=True)
            calc.optimize_orca(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.SMD,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            smd_orca_pos = calc.positions

            # 2 calculate Gibbs free energy with SMD
            calc.calculate_gibbs_free_energy(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.SMD,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            smd_free_energy = calc.gibbs_free_energy

            # 3 calculate Gibbs free energy with COSMORS
            calc.calculate_gibbs_free_energy(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.COSMORS,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            cosmo_rs_free_energy = calc.gibbs_free_energy

            # 3 calculate GNN free energy
            calc.calculate_gibbs_free_energy(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.CPCM,
                model=model,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            gnn_cpcm_free_energy_smd_opt = calc.gibbs_free_energy

            # Calculate Hessian of CPCM
            calc.calculate_hessian(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.CPCM,
                model=model,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            cpcm_hessian_at_smd_opt = calc.QM_hessian
            gnn_hessian_at_smd_opt = calc.GNN_hessian

            # Evaluate SMD plus GNN
            calc.calculate_gibbs_free_energy(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.SMD,
                model=model,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )
            gnn_smd_free_energy_smd_opt = calc.gibbs_free_energy

            # Evaluate single point energy
            calc.singlepoint(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.CPCM,
                model=model,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            cpcm_gnn_single_point_at_smd_opt = calc.electronic_energy

            # 4 Optimize with GNN
            calc.optimize(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.CPCM,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
                model=model,
                fmax=0.025,
            )

            calc.calculate_gibbs_free_energy(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.CPCM,
                model=model,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            gnn_opt_pos_0025 = calc.positions
            cpcm_gnn_free_energy_gnn_opt = calc.gibbs_free_energy

            # Get hessian at GNN optimized position
            calc.calculate_hessian(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.CPCM,
                model=model,
                nprocs=1,
                maxcore=10000,
                directory=create_temp_dir_with_random_name(base_tempdir),
            )

            cpcm_hessian_at_gnn_opt = calc.QM_hessian
            gnn_hessian_at_gnn_opt = calc.GNN_hessian

            # Evaluate single point energy
            calc.singlepoint(
                orca_solvent,
                Functional.B3LYP,
                BasisSet.def2_TZVP,
                SolventModel.CPCM,
                model=model,
                nprocs=1,
                maxcore=10000,
            )
            cpcm_gnn_single_point_at_gnn_opt = calc.electronic_energy

            attempts = 3
        except Exception as e:
            print(e)
            print("Retrying ...")
            print("copying directory to errors")
            os.makedirs(f"{args.output}_error/{args.solvent}_{args.idx}")
            os.system(
                f"cp -r {base_tempdir} {args.output}_error/{args.solvent}_{args.idx}"
            )
            attempts += 1

    result = (
        smd_orca_pos,
        smd_free_energy,
        cosmo_rs_free_energy,
        gnn_cpcm_free_energy_smd_opt,
        cpcm_hessian_at_smd_opt,
        gnn_hessian_at_smd_opt,
        gnn_smd_free_energy_smd_opt,
        cpcm_gnn_single_point_at_smd_opt,
        gnn_opt_pos_0025,
        cpcm_gnn_free_energy_gnn_opt,
        cpcm_hessian_at_gnn_opt,
        gnn_hessian_at_gnn_opt,
        cpcm_gnn_single_point_at_gnn_opt,
    )

    results.append(result)

np.save(
    f"{args.output}/{args.solvent}_{args.idx}.npy",
    np.array(
        results,
        dtype=object,
    ),
)
