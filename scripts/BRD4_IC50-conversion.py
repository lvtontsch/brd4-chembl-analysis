import pandas as pd
import numpy as np
from pathlib import Path

# extract
def extract_from_csv(source):
    dataframe = pd.read_csv(source, sep=';', low_memory=False)
    return dataframe

# transform
def transform(df):
    df['Standard Relation']=df['Standard Relation'].str.strip("'")
    
    # filter
    df = df[(df['Smiles'].notnull())
            & (df['Standard Type'] == 'IC50')
            & (df['Standard Relation'].isin(['=', '>']))
            & (df['Standard Value'].notnull())
            & (df['Standard Units'] == 'nM')
            & (df['Data Validity Comment'] != 'Potential transcription error')
            & (df['Potential Duplicate'] == 0)
            & (df['Assay Type'] == 'B')
            & (df['Target ChEMBL ID'] == 'CHEMBL1163125')
            & (df['Source Description'].isin(['BindingDB Patent Bioactivity Data', 'Scientific Literature']))
            ]
    keep = ['Molecule ChEMBL ID','Molecular Weight','#RO5 Violations','AlogP','Smiles','Standard Relation','Standard Value','pChEMBL Value','Assay ChEMBL ID','Assay Description','Assay Variant Mutation','Document ChEMBL ID','Source Description','Document Journal','Document Year']
    df = df[keep].copy()
    df['pIC50'] = (9 - np.log10(df['Standard Value'])).where(df['Standard Relation'] == '=')
    df['pIC50 upper bound'] = (9 - np.log10(df['Standard Value'])).where(df['Standard Relation'] == '>')

    # split mutations
    df_mutants = df[df['Assay Variant Mutation'].notna()]
    df_clean = df[df['Assay Variant Mutation'].isna()]
    
    # value per molecule
    df_molecules_exact = df_clean[df_clean['Standard Relation'] == '='].groupby('Molecule ChEMBL ID').agg(**{
        'Smiles': ('Smiles', 'first'),
        'Molecular Weight': ('Molecular Weight', 'first'),
        'AlogP': ('AlogP', 'first'),
        '#RO5 Violations': ('#RO5 Violations', 'first'),
        'pIC50': ('pIC50', 'median'),
        'pIC50 Spread': ('pIC50', lambda x: x.max() - x.min()), 
        'Number of Measurements': ('Standard Value', 'count'),
        'Number of Documents': ('Document ChEMBL ID', 'nunique'),
        'First Year': ('Document Year', 'min'),
        'Sources': ('Source Description', lambda x: ', '.join(sorted(x.unique())))
    }).reset_index()
    df_molecules_exact['Value Type'] = 'exact'
    df_molecules_exact['Activity'] = 'active'
    df_molecules_exact.loc[df_molecules_exact['pIC50'] < 6, 'Activity'] = 'inactive'
    df_molecules_exact['IC50 (nM)'] = 10 ** (9 - df_molecules_exact['pIC50'])

    df_molecules_lowerbound = df_clean[(df_clean['Standard Relation'] == '>') & (~df_clean['Molecule ChEMBL ID'].isin(df_molecules_exact['Molecule ChEMBL ID']))].groupby('Molecule ChEMBL ID').agg(**{
        'Smiles': ('Smiles', 'first'),
        'Molecular Weight': ('Molecular Weight', 'first'),
        'AlogP': ('AlogP', 'first'),
        '#RO5 Violations': ('#RO5 Violations', 'first'),
        'IC50 lower bound (nM)': ('Standard Value', 'max'),
        'Number of Measurements': ('Standard Value', 'count'),
        'Number of Documents': ('Document ChEMBL ID', 'nunique'),
        'First Year': ('Document Year', 'min'),
        'Sources': ('Source Description', lambda x: ', '.join(sorted(x.unique())))
    }).reset_index()
    df_molecules_lowerbound['Value Type'] = 'lower bound'
    df_molecules_lowerbound['Activity'] = 'unknown'
    df_molecules_lowerbound.loc[df_molecules_lowerbound['IC50 lower bound (nM)'] >= 1000, 'Activity'] = 'inactive'
    df_molecules_lowerbound = df_molecules_lowerbound[df_molecules_lowerbound['Activity'] == 'inactive'].copy()
    df_molecules_lowerbound['pIC50 upper bound'] = 9 - np.log10(df_molecules_lowerbound['IC50 lower bound (nM)'])

    df_molecules = pd.concat([df_molecules_exact, df_molecules_lowerbound], ignore_index=True)
    
    return df_clean, df_mutants, df_molecules

def load(df_clean, df_mutants, df_molecules, folder):
    folder.mkdir(parents=True, exist_ok=True)
    df_clean.to_csv(folder / "brd4_measurements.csv", index=False)
    df_mutants.to_csv(folder / "brd4_mutants.csv", index=False)
    df_molecules.to_csv(folder / "brd4_molecules.csv", index=False)

df = extract_from_csv(Path("data/raw/BRD4_IC50.csv"))
df_clean, df_mutants, df_molecules = transform(df)

load(df_clean, df_mutants, df_molecules, Path("data/processed/"))