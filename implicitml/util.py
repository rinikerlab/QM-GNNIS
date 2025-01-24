import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def embedMolecule(smiles):
    # TODO: check this function
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    return mol


def getAtomicNumbers(molecule):
    atomic_numbers = list()
    for atom in molecule.GetAtoms():
        atomic_numbers.append(atom.GetAtomicNum())
    return np.array(atomic_numbers, dtype=np.int64)


def normal_modes(hessian_matrix=None, realAtomIdxs=None, masses=None):
    """
    Adapted from https://leeping.github.io/forcebalance/doc/html/api/openmmio_8py_source.html

    OpenMM Normal Mode Analysis
    Since OpenMM doesnot provide a Hessian evaluation method, we used finite difference on forces

    Parameters
    ----------
    shot: int
        The frame number in the trajectory of this target
    optimize: bool, default True
        Optimize the geometry before evaluating the normal modes

    Returns
    -------
    freqs: np.array with shape (3N - 6) x 1, N = number of "real" atoms
        Harmonic frequencies, sorted from smallest to largest, with the 6 smallest removed, in unit cm^-1
    normal_modes: np.array with shape (3N - 6) x (3N), N = number of "real" atoms
        The normal modes corresponding to each of the frequencies, scaled by mass^-1/2.
    """

    # Mass weight hessian

    mass_weighted_hessian = np.zeros(hessian_matrix.shape)
    for i in range(len(realAtomIdxs) * 3):
        for j in range(len(realAtomIdxs) * 3):
            mass_weighted_hessian[i][j] = hessian_matrix[i][j] / np.sqrt(
                masses[i // 3] * masses[j // 3]
            )

    # step 1: build a full hessian matrix
    noa = len(realAtomIdxs)
    # step 2: diagonalize the hessian matrix
    eigvals, eigvecs = np.linalg.eigh(mass_weighted_hessian)
    # step 3: convert eigenvalues to frequencies
    coef = 0.5 / np.pi * 33.3564095  # 10^12 Hz => cm-1
    negatives = (eigvals >= 0).astype(int) * 2 - 1  # record the negative ones
    freqs = np.sqrt(eigvals + 0j) * coef * negatives
    # step 4: convert eigenvectors to normal modes
    # re-arange to row index and shape
    normal_modes = eigvecs.T.reshape(noa * 3, noa, 3)
    # step 5: Remove mass weighting from eigenvectors
    massList = np.array(masses)  # unit in dalton
    for i in range(normal_modes.shape[0]):
        mode = normal_modes[i]
        mode /= np.sqrt(massList[:, np.newaxis])
        mode /= np.linalg.norm(mode)
    # step 5: remove the 6 freqs with smallest abs value and corresponding normal modes
    n_remove = 5 if len(realAtomIdxs) == 2 else 6
    larger_freq_idxs = np.sort(np.argpartition(np.abs(freqs), n_remove)[n_remove:])
    # larger_freq_idxs = np.sort(np.argpartition(np.abs(freqs), n_remove))[n_remove:]
    freqs = freqs[larger_freq_idxs]
    normal_modes = normal_modes[larger_freq_idxs]
    return freqs, normal_modes
