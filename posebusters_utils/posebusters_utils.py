from posebusters import PoseBusters
from pathlib import Path
import streamlit as st
from rdkit import Chem
from utils.structures import Molecule

from posebusters.modules.intermolecular_distance import check_intermolecular_distance
from posebusters.tools.loading import safe_load_mol


def calculate_posebusters(mol_cond: Molecule, mol_true: Molecule, mol_pred: Molecule, mode='mol',**kwargs):
    # if not prediction_str:
    #     st.write('Please upload prediction sdf')
    #     return
    # if not protein_str and mode in ['redock', 'dock']:
    #     st.write('Please upload protein pdb')
    #     return
    # if not ground_truth_str and mode in ['redock']:
    #     st.write('Please upload ground truth sdf')
    #     return


    # print('suuup')
    # # TODO: don't save and load to file, transform string directly to rdkit mol
    # print(len(str_true), len(str_cond), len(str_pred))
    # mol_true = Chem.SDMolSupplier()
    # mol_true.SetData(str(str_true))
    # mol_true = mol_true[0]
    # print(mol_true)
    # mol_pred = Chem.SDMolSupplier()
    # mol_pred.SetData(str(str_pred))
    # mol_pred = mol_pred[0]
    # print(mol_pred)
    # mol_cond = Chem.rdmolfiles.MolFromPDBBlock(str(str_cond))
    # print(mol_true, mol_pred, mol_cond)


    # pb = PoseBusters(mode)
    # output = pb.bust(mol_pred=mol_pred, mol_true=mol_true, mol_cond=mol_cond)


    path_true = Path('tmp/ground_truth.sdf')
    path_pred = Path('tmp/prediction.sdf')
    path_cond = Path('tmp/protein.pdb')

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