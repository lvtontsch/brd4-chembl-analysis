# ChEMBL Analysis: BRD4 inhibitor viability

## Question
What makes a molecule a promising BRD4 inhibitor?

## Data
- Source: ChEMBL database (target: CHEMBL116125, IC50 measurements)
- Download date: 19.10.2026

## Data transformations
- Filtered to verify the sole inclusion of CHEMBL1163125 measurements of binding assays (type B) with IC50 values in nM
- Selection of exact and lower-bound measurements
- Exclusion of transcription errors and ChEMBL-flagged duplicates. 'Outside typical range' measurements were kept due to inconsistant application of the flag throughout the data.
- Selection of primary sources
- Extraction of relevant columns to minimise excessive data
- Calculation of pIC50 for measurements with exact results & calculation of 'pIC50 upper bound' for measurements with lower bound IC50 values
- Separation of mutant BRD4 measurements from the regular BRD4 data due to limited comparability within data of a selected protein
- Condensation of the cleaned data into two compartments of aggregates per molecule: exact measurements with pIC50 median calculation to properly respect the logarithmic relevance of the measured value & lower-bound measurements with max IC50 to 
- Classification into active and inactive inhibitors: threshold of pIC50 >= 6 (IC50 = 1 uM), common cutoff for active proteins

## Results
Dashboard: tba'

## Limitations
- 
- 

## How to run
1. Clone the repository and create the virtual environment:
    python -m venv venv
    .\venv\Scripts\activate
    pip install -r requirements.txt

2. Download the data from ChEMBL: select target CHEMBL1163125 (BRD4) and download the CSV of the IC50 activity assays. Save it as `data/raw/BRD4_IC50.csv`.

3. Run the data conversion script:
    python scripts/BRD4_IC50-conversoin.py

The output is written to `data/processed/`:
- `brd4_measurements.csv`: all cleaned measurements and calculation of pIC50
- `brd4_mutants.csv`: separated data on mutant BRD4
- `brd4_molecules.csv`: one row per molecule inclduing activity label