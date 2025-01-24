from enum import Enum, auto
from ase.calculators.orca import OrcaProfile, ORCA
from ase.calculators.mixing import SumCalculator
from ase.optimize import BFGS, LBFGS
from ase import Atoms
import shutil
from implicitml.util import getAtomicNumbers
from implicitml.ase import Calculator
from implicitml.units import (
    ENERGY_ASE_TO_IMPLICITML,
    FORCE_ASE_TO_IMPLICITML,
    LENGTH_ASE_TO_IMPLICITML,
    ENERGY_ORCA_TO_ASE,
    HESSIAN_ASE_TO_IMPLICITML,
)
import tempfile
from contextlib import nullcontext
from pathlib import Path
import logging
import pandas as pd
import torch
import numpy as np
import time

logger = logging.getLogger(__name__)


class Solvent(Enum):
    WATER = auto()
    CHLOROFORM = auto()
    BENZENE = auto()
    DMSO = auto()
    METHANOL = auto()
    THF = auto()
    ACETONITRILE = auto()
    PYRIDINE = auto()
    ACETONE = auto()
    ETHANOL = auto()
    ETHYLACETATE = auto()
    DCM = auto()
    DIETHYLETHER = auto()
    DIOXANE = auto()
    HEXANE = auto()

    def __str__(self):
        return self.name.lower()


class Functional(Enum):
    BP86 = auto()
    B3LYP = auto()
    PBE = auto()

    def __str__(self):
        return self.name


class BasisSet(Enum):
    def2_SVP = auto()
    def2_TZVP = auto()
    def2_TZVPD = auto()
    def2_QZVPP = auto()

    def __str__(self):
        return self.name.replace("_", "-")


class SolventModel(Enum):
    CPCM = auto()
    SMD = auto()
    COSMORS = auto()
    VACUUM = auto()

    def __str__(self):
        return self.name


DIELECTRIC = {
    Solvent.WATER: 78.5,
    Solvent.DMSO: 47.24,
    Solvent.CHLOROFORM: 4.81,
    Solvent.BENZENE: 2.28,
    Solvent.METHANOL: 33,
    Solvent.THF: 7.52,
    Solvent.ACETONITRILE: 36.64,
    Solvent.PYRIDINE: 13.26,
    Solvent.ACETONE: 21.01,
    Solvent.ETHANOL: 25.3,
    Solvent.ETHYLACETATE: 6.2,
    Solvent.DCM: 8.93,
    Solvent.DIETHYLETHER: 4.27,
    Solvent.DIOXANE: 2.22,
    Solvent.HEXANE: 1.88,
}


class ImplicitSolventCalculator:

    def __init__(self, molecule, orca_binary=None):
        if not orca_binary:
            orca_binary = shutil.which("orca")
            assert orca_binary is not None, "Provide a valid ORCA binary."
        assert (
            molecule.GetNumConformers() == 1
        ), "Molecule must be embedded with one conformer."
        symbols = getAtomicNumbers(molecule)
        positions = molecule.GetConformer(0).GetPositions()
        assert len(symbols) == len(positions)
        self.molecule = molecule
        self.atoms = Atoms(symbols, positions)
        self.profile = OrcaProfile(command=orca_binary)

    def singlepoint(
        self,
        solvent,
        functional,
        basis_set,
        solvent_model,
        dispersion_correction=True,
        model=None,
        directory=None,
        nprocs=1,
        maxcore=1024,
    ):
        assert solvent_model in [
            SolventModel.CPCM,
            SolventModel.SMD,
            SolventModel.VACUUM,
        ], f"Cannot do singlepoint calculation with {SolventModel.COSMORS}"
        with (
            tempfile.TemporaryDirectory() if not directory else nullcontext(directory)
        ) as directory:
            logger.info(f"ORCA singlepoint calculation in {directory}")
            inputs = [
                self._singlepoint_string(functional, basis_set, dispersion_correction),
            ]
            blocks = [
                self._resources_block(nprocs, maxcore),
            ]
            if solvent_model == SolventModel.CPCM:
                blocks.append(self._solvent_block(solvent, solvent_model))
            if solvent_model == SolventModel.SMD:
                inputs.append(self._solvent_string(solvent, solvent_model))
            calculators = list()
            orca_calc = ORCA(
                profile=self.profile,
                directory=directory,
                orcasimpleinput=" ".join(inputs),
                orcablocks="\n".join(blocks),
            )
            if model:
                torch_calc = Calculator(model)
                calculators.append(torch_calc)
            calculators.append(orca_calc)
            self.atoms.calc = SumCalculator(calculators)
            _ = self.atoms.get_potential_energy()

    def optimize_orca(
        self,
        solvent,
        functional,
        basis_set,
        solvent_model,
        dispersion_correction=True,
        model=None,
        directory=None,
        nprocs=1,
        maxcore=1024,
    ):
        assert solvent_model in [
            SolventModel.CPCM,
            SolventModel.SMD,
            SolventModel.VACUUM,
        ], f"Cannot do ORCA optimization calculation with {SolventModel.COSMORS}"
        assert model is None, "ORCA optimization does not support GNNs."
        with (
            tempfile.TemporaryDirectory() if not directory else nullcontext(directory)
        ) as directory:
            logger.info(f"ORCA Optimization calculation in {directory}")
            inputs = [
                self._orca_opt_string(functional, basis_set, dispersion_correction),
            ]
            blocks = [
                self._resources_block(nprocs, maxcore),
            ]
            if solvent_model == SolventModel.CPCM:
                blocks.append(self._solvent_block(solvent, solvent_model))
            if solvent_model == SolventModel.SMD:
                inputs.append(self._solvent_string(solvent, solvent_model))
            orca_calc = ORCA(
                profile=self.profile,
                directory=directory,
                orcasimpleinput=" ".join(inputs),
                orcablocks="\n".join(blocks),
            )
            self.atoms.calc = orca_calc
            _ = self.atoms.get_potential_energy()

            time.sleep(10)
            with open(f"{directory}/orca.out") as f:
                # check i "ORCA TERMINATED NORMALLY" in file
                if "ORCA TERMINATED NORMALLY" not in f.read():
                    print(f.read())
                    assert False, "ORCA calculation crashed."

            optimized_positions = pd.read_csv(
                f"{directory}/orca.xyz",
                delim_whitespace=True,
                skiprows=2,
                header=None,
                names=["atom", "x", "y", "z"],
            )

            self.atoms.positions = optimized_positions[["x", "y", "z"]].values

            # Check if set positions is available
            if hasattr(self.molecule.GetConformer(0), "SetPositions"):
                self.molecule.GetConformer(0).SetPositions(self.atoms.positions)
            else:
                from rdkit.Geometry import Point3D

                for i, pos in enumerate(self.atoms.positions):
                    conf = self.molecule.GetConformer(0)
                    conf.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))

    def optimize(
        self,
        solvent,
        functional,
        basis_set,
        solvent_model,
        dispersion_correction=True,
        model=None,
        directory=None,
        nprocs=1,
        maxcore=1024,
        fmax=0.05,
        maxstep=0.2,
        alpha=70.0,
        optimizer="BFGS",
    ):
        assert solvent_model in [
            SolventModel.CPCM,
            SolventModel.SMD,
            SolventModel.VACUUM,
        ], f"Cannot do geometry optmization with {SolventModel.COSMORS}"
        with (
            tempfile.TemporaryDirectory() if not directory else nullcontext(directory)
        ) as directory:
            logger.info(f"ORCA geometry optimization in {directory}")
            inputs = [
                self._optimization_string(functional, basis_set, dispersion_correction),
            ]
            blocks = [
                self._resources_block(nprocs, maxcore),
            ]
            if solvent_model == SolventModel.CPCM:
                blocks.append(self._solvent_block(solvent, solvent_model))
            if solvent_model == SolventModel.SMD:
                inputs.append(self._solvent_string(solvent, solvent_model))
            calculators = list()
            orca_calc = ORCA(
                profile=self.profile,
                directory=directory,
                orcasimpleinput=" ".join(inputs),
                orcablocks="\n".join(blocks),
            )
            if model:
                torch_calc = Calculator(model)
                calculators.append(torch_calc)
            calculators.append(orca_calc)
            self.atoms.calc = SumCalculator(calculators)
            if optimizer == "BFGS":
                opt = BFGS(self.atoms, maxstep=maxstep, alpha=alpha)
            elif optimizer == "LBFGS":
                opt = LBFGS(self.atoms, maxstep=maxstep, alpha=alpha)
            converged = opt.run(fmax=fmax)
            assert converged, "Geometry optimization has not converged!"

            # Check if set positions is available
            if hasattr(self.molecule.GetConformer(0), "SetPositions"):
                self.molecule.GetConformer(0).SetPositions(self.atoms.positions)
            else:
                from rdkit.Geometry import Point3D

                for i, pos in enumerate(self.atoms.positions):
                    conf = self.molecule.GetConformer(0)
                    conf.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))

    def calculate_hessian(
        self,
        solvent,
        functional,
        basis_set,
        solvent_model,
        dispersion_correction=True,
        directory=None,
        model=None,
        nprocs=1,
        maxcore=1024,
    ):

        assert solvent_model in [
            SolventModel.CPCM,
            SolventModel.SMD,
            SolventModel.VACUUM,
        ], f"Cannot calculate Hessian for {SolventModel.COSMORS}"

        with (
            tempfile.TemporaryDirectory() if not directory else nullcontext(directory)
        ) as directory:
            self.calculate_gibbs_free_energy(
                solvent,
                functional,
                basis_set,
                solvent_model,
                dispersion_correction,
                directory,
                model,
                nprocs,
                maxcore,
            )

            time.sleep(10)
            hessian = self._parse_hessian(Path(directory) / "orca.hess")
            self.atoms.calc.results["QM_hessian"] = np.array(hessian)

        if model:
            positions = self.positions
            tensor_positions = torch.tensor(
                positions, dtype=torch.float32, requires_grad=True
            )
            hessian = torch.autograd.functional.hessian(model, tensor_positions)
            self.atoms.calc.results["GNN_hessian"] = np.array(
                hessian.view(hessian.shape[0] * 3, hessian.shape[0] * 3)
                / HESSIAN_ASE_TO_IMPLICITML
            )

    def calculate_gibbs_free_energy(
        self,
        solvent,
        functional,
        basis_set,
        solvent_model,
        dispersion_correction=True,
        directory=None,
        model=None,
        nprocs=1,
        maxcore=1024,
    ):
        with (
            tempfile.TemporaryDirectory() if not directory else nullcontext(directory)
        ) as directory:
            if solvent_model in [
                SolventModel.CPCM,
                SolventModel.SMD,
                SolventModel.VACUUM,
            ]:
                logger.info(f"ORCA Gibbs free energy calculation in {directory}")
                inputs = [
                    self._free_energy_string(
                        functional, basis_set, dispersion_correction
                    ),
                ]
                blocks = [
                    self._resources_block(nprocs, maxcore),
                ]
                if solvent_model == SolventModel.CPCM:
                    blocks.append(self._solvent_block(solvent, solvent_model))
                if solvent_model == SolventModel.SMD:
                    inputs.append(self._solvent_string(solvent, solvent_model))
                calculators = list()
                orca_calc = ORCA(
                    profile=self.profile,
                    directory=directory,
                    orcasimpleinput=" ".join(inputs),
                    orcablocks="\n".join(blocks),
                )
                if model:
                    torch_calc = Calculator(model)
                    calculators.append(torch_calc)
                calculators.append(orca_calc)
                self.atoms.calc = SumCalculator(calculators)
                _ = self.atoms.get_potential_energy()
                gibbs_free_energy = self._parse_gibbs_free_energy(
                    Path(directory) / "orca.out"
                )
                gibbs_free_energy *= ENERGY_ORCA_TO_ASE
                self.atoms.calc.results["free_energy"] = gibbs_free_energy
                if model:
                    self.atoms.calc.results["free_energy"] += torch_calc.results[
                        "energy"
                    ]

            elif solvent_model == SolventModel.COSMORS:
                logger.info(f"ORCA Gibbs free energy calculation in {directory}")
                # vacuum calculation, without GNN
                inputs = [
                    self._free_energy_string(
                        functional, basis_set, dispersion_correction
                    ),
                ]
                blocks = [
                    self._resources_block(nprocs, maxcore),
                ]
                orca_calc = ORCA(
                    profile=self.profile,
                    directory=directory,
                    orcasimpleinput=" ".join(inputs),
                    orcablocks="\n".join(blocks),
                )
                self.atoms.calc = orca_calc
                _ = self.atoms.get_potential_energy()
                gibbs_free_energy = self._parse_gibbs_free_energy(
                    Path(directory) / "orca.out"
                )

                # COSMO-RS calculation, with GNN if requested (would mean GNN on COSMO-RS!)
                inputs = [
                    self._solvent_string(solvent, solvent_model),
                ]
                blocks = [
                    self._resources_block(nprocs, maxcore),
                ]
                calculators = list()
                orca_calc = ORCA(
                    profile=self.profile,
                    directory=directory,
                    orcasimpleinput=" ".join(inputs),
                    orcablocks="\n".join(blocks),
                )
                if model:
                    torch_calc = Calculator(model)
                    calculators.append(torch_calc)
                calculators.append(orca_calc)
                self.atoms.calc = SumCalculator(calculators)
                _ = self.atoms.get_potential_energy()
                solvation_free_energy = self._parse_solvation_free_energy(
                    Path(directory) / "orca.out"
                )

                # Gsolv = Gvac + dGsolv
                self.atoms.calc.results["free_energy"] = (
                    gibbs_free_energy + solvation_free_energy
                ) * ENERGY_ORCA_TO_ASE
            else:
                raise ValueError(f"Unknown solvent model: {solvent_model}")

    @property
    def electronic_energy(self):
        electronic_energy = self.atoms.calc.results.get("energy")
        if electronic_energy:
            electronic_energy *= ENERGY_ASE_TO_IMPLICITML
        else:
            logger.warning(f"No entry found for energy. Have you calculated it?")
        return electronic_energy

    @property
    def gibbs_free_energy(self):
        gibbs_free_energy = self.atoms.calc.results.get("free_energy")
        if gibbs_free_energy:
            gibbs_free_energy *= ENERGY_ASE_TO_IMPLICITML
            if gibbs_free_energy == self.electronic_energy:
                logger.warning(
                    f"Electronic energy and Gibbs free energy are identical."
                )
        else:
            logger.warning(f"No entry found for free energy. Have you calculated it?")
        return gibbs_free_energy

    @property
    def forces(self):
        forces = self.atoms.calc.results.get("forces")
        if forces is not None:
            forces *= FORCE_ASE_TO_IMPLICITML
        else:
            logger.warning(f"No entry found for forces. Have you calculated them?")
        return forces

    @property
    def positions(self):
        assert (
            self.molecule.GetConformer(0).GetPositions() == self.atoms.positions
        ).all()
        return self.molecule.GetConformer(0).GetPositions() * LENGTH_ASE_TO_IMPLICITML

    @property
    def hessian(self):

        hessian = None
        if self.atoms.calc.results.get("QM_hessian") is not None:
            hessian = self.QM_hessian
        if self.atoms.calc.results.get("GNN_hessian") is not None:
            if hessian is not None:
                hessian += self.GNN_hessian
            else:
                hessian = self.GNN_hessian
        if hessian is None:
            logger.warning(f"No entry found for hessian. Have you calculated it?")
        return hessian

    @property
    def QM_hessian(self):
        hessian = self.atoms.calc.results.get("QM_hessian")
        if not hessian is None:
            hessian = np.array(hessian) * HESSIAN_ASE_TO_IMPLICITML
        else:
            logger.warning(f"No entry found for hessian. Have you calculated it?")
        return hessian

    @property
    def GNN_hessian(self):
        hessian = self.atoms.calc.results.get("GNN_hessian")
        if not hessian is None:
            hessian = np.array(hessian) * HESSIAN_ASE_TO_IMPLICITML
        else:
            logger.warning(f"No entry found for hessian. Have you calculated it?")
        return hessian

    def _singlepoint_string(self, functional, basis_set, dispersion_correction):
        assert functional in Functional, f"Invalid functional: {functional}"
        assert basis_set in BasisSet, f"Invalid basis set: {basis_set}"
        line = f"{functional} {basis_set}"
        if dispersion_correction:
            line += " D4"
        return line

    def _orca_opt_string(self, functional, basis_set, dispersion_correction):
        assert functional in Functional, f"Invalid functional: {functional}"
        assert basis_set in BasisSet, f"Invalid basis set: {basis_set}"
        line = f"{functional} {basis_set} OPT"
        if dispersion_correction:
            line += " D4"
        return line

    def _optimization_string(self, functional, basis_set, dispersion_correction):
        assert functional in Functional, f"Invalid functional: {functional}"
        assert basis_set in BasisSet, f"Invalid basis set: {basis_set}"
        singlepoint_string = self._singlepoint_string(
            functional, basis_set, dispersion_correction
        )
        return singlepoint_string + " " + "EnGrad"

    def _free_energy_string(self, functional, basis_set, dispersion_correction):
        assert functional in Functional, f"Invalid functional: {functional}"
        assert basis_set in BasisSet, f"Invalid basis set: {basis_set}"
        singlepoint_string = self._singlepoint_string(
            functional, basis_set, dispersion_correction
        )
        return singlepoint_string + " " + "Freq"

    def _solvent_string(self, solvent, solvent_model):
        assert solvent in Solvent, f"Invalid solvent: {solvent}"
        assert solvent_model in SolventModel, f"Invalid solvent model: {solvent_model}"
        assert (
            solvent_model != SolventModel.CPCM
        ), "Use solvent block for CPCM to include epsilon."
        return f"{solvent_model}({solvent})"

    def _resources_block(self, nprocs, maxcore):
        return f"! MINIPRINT \n %pal nprocs {nprocs} end %maxcore {maxcore}"

    def _solvent_block(self, solvent, solvent_model):
        assert (
            solvent_model == SolventModel.CPCM
        ), "Use solvent string for any model that is not CPCM."
        return f"%cpcm epsilon {DIELECTRIC[solvent]} end"

    def _parse_gibbs_free_energy(self, file):
        with open(file) as f:
            for line in f:
                if "Final Gibbs free energy" in line:
                    return float(line.split()[5])

    def _parse_solvation_free_energy(self, file):
        with open(file) as f:
            for line in f:
                if "Free energy of solvation" in line:
                    return float(line.split()[6])

    def _parse_hessian(self, file):

        file_lines = open(file, "r", encoding="utf-8").readlines()

        hessian_blocks = []

        for i, line in enumerate(file_lines):
            if "$hessian" not in line:
                continue

            # Ensure the number of atoms is present, and is the number expected
            n_atoms = int(file_lines[i + 1].split()[0]) // 3
            break

        start_line = i + 3

        for j, h_line in enumerate(file_lines[start_line:]):
            if len(h_line.split()) == 0:
                # Assume we're at the end of the Hessian
                break

            # Skip blank lines in the file, marked by one or more fewer items
            # than the previous
            if len(h_line.split()) < len(file_lines[start_line + j - 1].split()):
                continue

            # First item is the coordinate number, thus append all others
            hessian_blocks.append([float(v) for v in h_line.split()[1:]])
            if "#" in h_line:
                break

        n_atoms = self.positions.shape[0]
        hessian = [block for block in hessian_blocks[: 3 * n_atoms]]

        for i, block in enumerate(hessian_blocks[3 * n_atoms :]):
            hessian[i % (3 * n_atoms)] += block

        return hessian
