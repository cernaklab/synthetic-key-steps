from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Atom, BondType

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def amat_entry(amat,r,c,val):
    """
    changes the (r,c) and (c,r) values of matrix "amat" to "val"
    """
    amat[r-1][c-1] = amat[c-1][r-1]= val
    return


def amat_edit(amat,r,c,delta):
    """
    changes the (r,c) and (c,r) values of matrix "amat" to "val"
    """
    amat[r-1][c-1] += delta
    if r != c:
        amat[c-1][r-1] += delta
    return


def molFromAdjMat(atoms, amat,sanitize=True):
    """Creates a mol object from an adjacency matrix.
    Inputs:
    atoms: list of atomic numbers of atoms, by row
    amat: adjacency matrix. Has to have same length as atoms (obviously)
    Output: mol object
    """
    
    m = Chem.RWMol()
    # add in the separate atoms
    for a in atoms: m.AddAtom(Atom(a))
    side_len = len(amat)    
    for r in range(side_len):
        for c in range(r+1,side_len):
            bond_order = amat[r][c]
            if bond_order > 0:
                if bond_order == 1: m.AddBond(r,c,BondType.SINGLE)
                if bond_order == 2: m.AddBond(r,c,BondType.DOUBLE)
                if bond_order == 3: m.AddBond(r,c,BondType.TRIPLE)

    if sanitize:
        Chem.SanitizeMol(m)
    return m


def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx()+1)
    return mol

def make_changelogs(data_file_path):
    data = pd.read_csv(data_file_path)[["bond","edit","file"]]
    data = data[~data.bond.isnull()].copy()
    
    changelogs = []
    entry = {}
    bond_edits = ()
    for r in data.itertuples():


        # reset entry dict at new step
        if r[1] == "step":
            entry["edits"] = bond_edits
            changelogs.append(entry)
            entry = {}

        # if padding, fill th
        elif r[1] == "pad":
            pad_atoms = r[3].split(" ")
            pad_atoms = [int(i) for i in pad_atoms]

            if pad_atoms == [0]:
                entry["pad"] = 0
                entry["pad_elem"] = []

            else:
                entry["pad"] = len(pad_atoms)
                entry["pad_elem"] = pad_atoms

            # make the empty bond edits here to prepare
            bond_edits = []
        else:
            bond_edits.append((int(r[1]),int(r[2]),int(r[3])))
    
    return changelogs


def apply_changes(amat_init, atoms_init,changelogs):
    seq_out = [amat_init.copy()]
    amat = amat_init.copy()
    atoms = atoms_init.copy()
    for i in changelogs:
        try:
            
            pad_amt = i["pad"]

            if pad_amt > 0:
                amat = np.pad(amat,[(0, pad_amt), (0, pad_amt)],  mode="constant")
                atoms.extend(i["pad_elem"])

            for ed in i["edits"]:
                amat_edit(amat,ed[0],ed[1],ed[2])
            seq_out.append(amat.copy())
            
        except:
            print(i)
        
    seq_out.reverse()
    
    all_sizes = [m.shape[0] for m in seq_out]
    max_size = max(all_sizes)

    output_padded = []

    for mat in seq_out:
        mat_size = mat.shape[0]
        if mat_size < max_size:
            pad_size = max_size - mat_size 
            output_padded.append(np.pad(mat, [(0, pad_size), (0, pad_size)], mode='constant'))
        else:
            output_padded.append(mat)
            
        
    return output_padded,atoms
    