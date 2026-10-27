# ChEMBL Analysis: BRD4 inhibitor potency

## Question
BRD4 is an epigenetic reader protein that drives oncogene expression (e.g. MYC), making it a target for cancer therapy.
How drug-like are potent BRD4 inhibitors, and how has the discovery of new inhibitors evolved since 2011?

## Data
- Source: ChEMBL (EMBL-EBI); target: CHEMBL1163125, IC50 measurements
- Download date: 19.10.2026

## Data transformations
- Filtered to verify the sole inclusion of CHEMBL1163125 measurements of binding assays (type B) with IC50 values in nM
- Selection of exact (=) measurements only; < and > values were excluded as they provide no exact potency value; < values originated almost entirely from patent potency categories (e.g. < 100 nM)
- Exclusion of transcription errors and ChEMBL-flagged duplicates. 'Outside typical range' measurements were kept due to inconsistent application of the flag throughout the data
- Exclusion of (−)-JQ1 (CHEMBL3220211), the inactive enantiomer used as a negative control, due to contradictory values (likely mis-assigned stereochemistry)
- Selection of primary sources (scientific literature and patent data); secondary compilations (EUbOPEN, SGC, chemical probe sets) were excluded, as they partly re-report previously published values and contributed < 60 measurements
- Extraction of relevant columns to minimise excessive data
- Reclassification of assay format (derived from ChEMBL's BAO label, corrected for protein expression hosts Sf9, BL21 and CHO) and annotation of assay technology (Alpha, TR-FRET, FP) 
- Calculation of pIC50
- Separation of mutant BRD4 measurements from the regular BRD4 data as mutations may alter the binding site and potencies are not directly comparable to wild-type BRD4
- Condensation of protein assay-based measurements into aggregates per molecule: cell-, tissue- and organism-based values were excluded due to permeability and downstream effects on IC50 results; key information is the median pIC50, aggregation on the logarithmic scale weights fold-differences equally and the median limits the influence of outliers
- Classification into active and inactive inhibitors: threshold of pIC50 >= 6 (IC50 <= 1 µM), a common cutoff for active compounds

## Results
[![BRD4 dashboard](https://public.tableau.com/static/images/br/brd4_analysis/Dashboard/1.png)](https://public.tableau.com/views/brd4_analysis/Dashboard)

Interactive version: [Tableau Public](https://public.tableau.com/views/brd4_analysis/Dashboard)

Molecules without Rule-of-5 data in ChEMBL (n = 25) are not shown in the Rule-of-5 charts.

## Limitations
- Measurements on different BRD4 constructs (BD1, BD2, tandem) were combined, as the domain could not be reliably determined from assay descriptions; for domain-selective compounds, the median pIC50 therefore represents neither domain (the largest spreads in the data result from BD2-selective inhibitors)
- Differences in assay conditions between laboratories may be reflected in the spread of repeat measurements of the same molecule
- Actives outnumber inactives (4,107 vs. 2,934), likely reflecting publication bias towards optimised compound series
- Accumulation of pIC50 values at 4.3 (50 µM) suggests assay upper limits recorded as exact values; these were retained, as they are classified as inactive
- Stereochemical annotation may be unreliable, so that measurements of one enantiomer may be assigned to the other

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
- `brd4_measurements.csv`: all cleaned exact measurements and calculation of pIC50
- `brd4_mutants.csv`: separated data on mutant BRD4
- `brd4_molecules.csv`: one row per molecule including activity label
