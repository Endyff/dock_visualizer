from posebusters import PoseBusters
from pathlib import Path
import streamlit as st
from rdkit import Chem
from utils.structures import Molecule

from posebusters.modules.intermolecular_distance import check_intermolecular_distance
from posebusters.tools.loading import safe_load_mol


def calculate_posebusters(mol_cond: Molecule, mol_true: Molecule, mol_pred: Molecule, mode='mol',**kwargs):
   
    tmp = Path('tmp')
    tmp.mkdir(exist_ok=True)
    path_true = tmp / Path('ground_truth.sdf')
    path_pred = tmp / Path('prediction.sdf')
    path_cond = tmp / Path('protein.pdb')

    with open(path_true, 'w') as f:
        f.write(str(mol_true))
    with open(path_pred, 'w') as f:
        f.write(str(mol_pred))
    with open(path_cond, 'w') as f:
        f.write(str(mol_cond))

    pb = PoseBusters(mode)
    output = pb.bust(mol_pred=path_pred, mol_true=path_true, mol_cond=path_cond)
    # print(output)
    fail_checks = []
    # print(output.T)
    for c in output.columns:
        try:
            if output[c].values[0] == False:
                fail_checks.append(f'{c}\n')
        except Exception as e:
            print(e)
            pass
    return fail_checks

def check_distance(str_cond, str_pred, radius_type):

    # mol_true = Path('tmp/ground_truth.sdf')
    mol_pred = Path('tmp/prediction.sdf')
    mol_cond = Path('tmp/protein.pdb')

    # with open(mol_true, 'w') as f:
    #     f.write(str(str_true))
    with open(mol_pred, 'w') as f:
        f.write(str(str_pred))
    with open(mol_cond, 'w') as f:
        f.write(str(str_cond))
    
    str_cond = safe_load_mol(mol_cond)
    # str_true = safe_load_mol(mol_true)
    str_pred = safe_load_mol(mol_pred)
    return check_intermolecular_distance(str_pred, str_cond, radius_type)