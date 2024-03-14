import py3Dmol
import pandas as pd
import streamlit as st
import argparse
from stmol import showmol

from utils.pdb_utils import PDBMolecule
from utils.sdf_utils import SDFMolecule
from posebusters_utils.posebusters_utils import calculate_posebusters, check_distance
from mol_visualizer import *
from pathlib import Path

from collections import Counter


parser = argparse.ArgumentParser()
parser.add_argument('--csv', help='Path to CSV file with paths')
parser.add_argument('--height', help='Height of the 3Dmol viewer', default=600)
parser.add_argument('--width', help='Width of the 3Dmol viewer', default=800)
args = parser.parse_args()

# Access the CSV file path using args.csv
df = None
if args.csv:
    df = pd.read_csv(args.csv)
else:
    df = pd.read_csv(st.sidebar.file_uploader('Upload CSV', type='csv'))
if df is not None:
    view = py3Dmol.view(width=args.width, height=args.height)
    complex = st.sidebar.selectbox('Select complex', df['complex_name'])
    row = df.loc[df['complex_name'] == complex]

    with open(row['mol_pred'].iloc[0]) as f:
        mol_pred = SDFMolecule(f.read())
    if 'mol_cond' in row and row['mol_cond'].iloc[0]:
        if Path(row['mol_cond'].iloc[0]).exists():
            with open(row['mol_cond'].iloc[0]) as f:
                mol_cond = PDBMolecule(f.read())
        else:
            st.write(f'File {row["mol_cond"].iloc[0]} does not exist')

    if 'mol_true' in row and row['mol_true'].iloc[0]:
        if Path(row['mol_true'].iloc[0]).exists():
            with open(row['mol_true'].iloc[0]) as f:
                mol_true = SDFMolecule(f.read())
        else:
            st.write(f'File {row["mol_true"].iloc[0]} does not exist')


    # RMSD        
    center_true = mol_true.get_mean()
    center_pred = mol_pred.get_mean()
    center_rmsd = ((center_true['x']-center_pred['x'])**2 + (center_true['y']-center_pred['y'])**2 + (center_true['z']-center_pred['z'])**2)**0.5
    st.write(f'Center of mass RMSD: {center_rmsd:.2f} Å')

    rmsd = mol_pred.rmsd(mol_true)
    if type(rmsd) == float:
        st.write(f'RMSD: {center_rmsd:.2f} Å')
    else:
        c1 = Counter([a['element'] for a in mol_pred.atoms])
        c2 = Counter([a['element'] for a in mol_pred.atoms])
        err_msg = " Check if the molecules have the same number of atoms" if c1 != c2 else ""
        st.write(f'RMSD calculation failed.{err_msg}')
        
        print(f'RMSD calculation failed, error message: {rmsd}')

    kwargs = {}

    hide_cond = st.sidebar.checkbox('Hide target', value=False)
    hide_true = st.sidebar.checkbox('Hide ground truth (green)', value=False)
    hide_pred = st.sidebar.checkbox('Hide prediction (blue)', value =False)


    pb_mode = st.sidebar.selectbox('Posebusters mode', ["redock", "dock", "mol"])
    pb_output = calculate_posebusters(mol_cond=mol_cond,
                                    mol_true=mol_true,
                                    mol_pred=mol_pred,
                                    mode=pb_mode)



    show_clash_results = st.sidebar.checkbox('Show clash results', value=True)
    show_clash_radii = st.sidebar.checkbox('Show clash radii', value=False)
    radius_type = st.sidebar.selectbox('Clash radii type', ['vdw', 'covalent'])

    show_weak_clashes = st.sidebar.checkbox('Show weak clashes', value=False)
    if show_clash_results:
        result_df = check_distance(mol_cond, mol_pred, radius_type)
        result_df['details'] = result_df['details'].loc[result_df['details'].clash]
        for i, row in result_df['details'].iterrows():
            at0 = mol_pred[row['ligand_atom_id']]
            at1 = mol_cond[row['protein_atom_id'] + 1]

            if row["protein_element"] == 'H':  # ignore clashes with hydrogens
                continue
            res_ats = mol_cond.select_residue(at1)
            mol_cond.set_residue_style(at1, {'stick':{'radius':0.1}})
            draw_line(view, at0, at1, {'color':'red'})

            if show_clash_radii:
                draw_sphere(view, at1, {'color':'red', 'radius': row['protein_vdw'], 'opacity': 0.7})
                draw_sphere(view, at0, {'color':'blue', 'radius': row['ligand_vdw'], 'opacity': 0.5}) 

    show_lines = st.sidebar.checkbox('Show lines', value=False)

    if show_lines:
        kwargs['line'] = {'line': {}}

    create_view(view, 
                mol_true=mol_true,
                mol_cond=mol_cond,
                mol_pred=mol_pred,
                hide_true=hide_true,
                hide_cond=hide_cond,
                hide_pred=hide_pred,
                show_lines=show_lines,
                )

    st.sidebar.write('PoseBusters failed checks:')
    for line in pb_output:
        st.sidebar.write(line)
    showmol(view, width=args.width, height=args.height)