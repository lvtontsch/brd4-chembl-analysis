# ChEMBL Analysis: BRD4 inhibitor potency

## Question
What makes a molecule a promising BRD4 inhibitor?

## Data
- Source: ChEMBL database (target: CHEMBL1163125, IC50 measurements)
- Download date: 19.10.2026

## Data transformations
- Filtered to verify the sole inclusion of CHEMBL1163125 measurements of binding assays (type B) with IC50 values in nM
- Selection of exact (=) and lower-bound (>) measurements; < measurements (~3,900) were excluded, as they provide no exact potency value and originated almost entirely from patents reporting potency categories (e.g. < 100 nM) rather than measured values
- Exclusion of transcription errors and ChEMBL-flagged duplicates. 'Outside typical range' measurements were kept due to inconsistent application of the flag throughout the data
- Selection of primary sources (scientific literature and patent data); secondary compilations (EUbOPEN, SGC, chemical probe sets) were excluded, as they partly re-report previously published values and contributed < 60 measurements
- Extraction of relevant columns to minimise excessive data
- Calculation of pIC50 for measurements with exact results & calculation of 'pIC50 upper bound' for measurements with lower bound IC50 values
- Separation of mutant BRD4 measurements from the regular BRD4 data as mutations may alter the binding site and potencies are not directly comparable to wild-type BRD4
- Condensation of the cleaned data into two compartments of aggregates per molecule: exact measurements via median pIC50, as aggregation on the logarithmic scale weights equal fold-differences equally and the median limits the influence of outliers; lower-bound measurements via maximum IC50, retaining the strongest constraint
- Classification into active and inactive inhibitors: threshold of pIC50 ≥ 6 (IC50 ≤ 1 µM), a common cutoff for active compounds; lower-bound molecules with IC50 > 1 µM were classified as inactive, those with lower bounds < 1 µM were excluded as uninformative

## Results
Dashboard: tba

## Limitations
- Measurements on different BRD4 constructs (BD1, BD2, tandem) were combined, as the domain could not be reliably determined from assay descriptions
- 358 molecules show a pIC50 spread > 1 (> 10-fold difference) between repeat measurements
- Lower-bound (`>`) values were used as inactive labels, while `<` values were excluded
- Actives outnumber inactives (4,449 vs. 3,606), likely reflecting publication bias towards optimised compound series

## How to run
1. Clone the repository and create the virtual environment:
    ```
    python -m venv venv
    .\venv\Scripts\activate 
    pip install -r requirements.txt
    ```

2. Download the data from ChEMBL: select target CHEMBL1163125 (BRD4) and download the CSV of the IC50 activity assays. Save it as `data/raw/BRD4_IC50.csv`.

3. Run the data conversion script:
    ```
    python scripts/BRD4_IC50-conversion.py
    ```

The output is written to `data/processed/`:
- `brd4_measurements.csv`: all cleaned measurements and calculation of pIC50
- `brd4_mutants.csv`: separated data on mutant BRD4
- `brd4_molecules.csv`: one row per molecule including activity label