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
            & (df['Standard Relation'] == '=')
            & (df['Standard Value'].notnull())
            & (df['Standard Units'] == 'nM')
            & (df['Data Validity Comment'] != 'Potential transcription error')
            & (df['Potential Duplicate'] == 0)
            & (df['Assay Type'] == 'B')
            & (df['Target ChEMBL ID'] == 'CHEMBL1163125')
            & (df['Source Description'].isin(['BindingDB Patent Bioactivity Data', 'Scientific Literature']))
            ]
    
    # select columns
    keep = ['Molecule ChEMBL ID','Molecular Weight','#RO5 Violations','AlogP','Smiles','Standard Relation','Standard Value','pChEMBL Value','Assay ChEMBL ID','Assay Description','BAO Label','Assay Cell Type','Assay Variant Mutation','Document ChEMBL ID','Source Description','Document Journal','Document Year']
    df = df[keep]

    # calculate pIC50 for exact Standard Values
    df['pIC50'] = (9 - np.log10(df['Standard Value']))
    
    # add assay categories
    df['Assay Technology'] = 'other/unclear'
    df.loc[df['Assay Description'].str.contains('alphascreen|alphalisa|alpha screen|alpha-screen', case=False, na=False), 'Assay Technology'] = 'Alpha'
    df.loc[df['Assay Description'].str.contains('FRET|HTRF|europium|Eu-labeled|terbium|allophycocyanin', case=False, na=False), 'Assay Technology'] = 'TR-FRET'
    df.loc[df['Assay Description'].str.contains('polari', case=False, na=False), 'Assay Technology'] = 'FP'
    df['Assay Format'] = df['BAO Label']
    df.loc[df['Assay Cell Type'].isin(['Sf9', 'BL21-DE3', 'CHO']), 'Assay Format'] = 'single protein format'

    # split mutations
    df_mutants = df[df['Assay Variant Mutation'].notna()].copy()
    df_clean = df[df['Assay Variant Mutation'].isna()].copy()
    
    # value per molecule
    df_molecules = df_clean[(df_clean['Assay Format'].isin(['single protein format', 'assay format']))].groupby('Molecule ChEMBL ID').agg(**{
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
    df_molecules['Activity'] = 'active'
    df_molecules.loc[df_molecules['pIC50'] < 6, 'Activity'] = 'inactive'
    df_molecules['IC50 (nM)'] = 10 ** (9 - df_molecules['pIC50'])
    
    return df_clean, df_mutants, df_molecules

def load(df_clean, df_mutants, df_molecules, folder):
    folder.mkdir(parents=True, exist_ok=True)
    df_clean.to_csv(folder / "brd4_measurements.csv", index=False)
    df_mutants.to_csv(folder / "brd4_mutants.csv", index=False)
    df_molecules.to_csv(folder / "brd4_molecules.csv", index=False)

df = extract_from_csv(Path("data/raw/BRD4_IC50.csv"))
df_clean, df_mutants, df_molecules = transform(df)

load(df_clean, df_mutants, df_molecules, Path("data/processed/"))