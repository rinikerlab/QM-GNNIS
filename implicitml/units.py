from ase import units as u

# (GNN -> ASE)
ENERGY_TORCH_TO_ASE = u.kJ / u.mol
FORCE_TORCH_TO_ASE = u.kJ / (u.mol * u.nm)

# (ASE -> GNN)
LENGTH_ASE_TO_TORCH = 1 / u.nm

# (ORCA output -> ASE)
ENERGY_ORCA_TO_ASE = u.Hartree

# (ASE -> implicitml calculator)
ENERGY_ASE_TO_IMPLICITML = 1 / (u.kJ / u.mol)
FORCE_ASE_TO_IMPLICITML = 1 / (u.kJ / (u.mol * u.nm))
HESSIAN_ASE_TO_IMPLICITML = 1 / (u.kJ / (u.mol * u.nm**2))
LENGTH_ASE_TO_IMPLICITML = 1 / (u.nm)
