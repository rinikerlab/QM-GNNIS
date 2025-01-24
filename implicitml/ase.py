# -*- coding: utf-8 -*-
"""Tools for interfacing with `ASE`_.

.. _ASE:
    https://wiki.fysik.dtu.dk/ase
"""

import torch
import ase.calculators.calculator
from implicitml.units import (
    ENERGY_TORCH_TO_ASE,
    FORCE_TORCH_TO_ASE,
    LENGTH_ASE_TO_TORCH,
)


class Calculator(ase.calculators.calculator.Calculator):

    implemented_properties = ["energy", "forces"]

    def __init__(self, model):
        super().__init__()
        self.model = model
        # no gradients on model parameters are required here
        for p in self.model.parameters():
            p.requires_grad_(False)

        a_parameter = next(self.model.parameters())
        self.device = a_parameter.device
        self.dtype = a_parameter.dtype

    def calculate(
        self,
        atoms=None,
        properties=["energy"],
        system_changes=ase.calculators.calculator.all_changes,
    ):
        super().calculate(atoms, properties, system_changes)
        print("Calculating TORCH")
        coordinates = torch.tensor(self.atoms.get_positions())
        coordinates = (
            coordinates.to(self.device)
            .to(self.dtype)
            .requires_grad_("forces" in properties)
        )
        coordinates = coordinates * LENGTH_ASE_TO_TORCH
        energy = self.model(coordinates)

        # Forces need to be calculated before energy is converted to eV otherwise the gradient will be wrong!
        if "forces" in properties:
            forces = -torch.autograd.grad(energy.squeeze(), coordinates)[0]
            forces = forces * FORCE_TORCH_TO_ASE
            self.results["forces"] = forces.squeeze(0).to("cpu").numpy()

        # native units of ASE are eV, Angstrom, and K
        energy = energy * ENERGY_TORCH_TO_ASE
        self.results["energy"] = energy.item()
