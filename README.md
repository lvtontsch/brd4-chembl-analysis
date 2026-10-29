# ChEMBL Analysis: BRD4 inhibitor potency
BRD4 is an epigenetic reader protein that drives oncogene expression (e.g. MYC), making it a target for cancer therapy.

## Question
How drug-like are potent BRD4 inhibitors and how have the number and drug-likeness of new inhibitors evolved since 2011?

## Data
- Source: ChEMBL 37 (EMBL-EBI), SQLite database
- Target: BRD4 (CHEMBL1163125), IC50 measurements
- Download date: 29.10.2026

## Data transformations
### Extraction and filtering: SQL, `sql/extract_brd4.sql`
Queries supporting the filtering decisions: `notebooks/filter_checks.ipynb`

- Joining of activity, assay, molecule, compound property, document and source tables
- Filtered to verify the sole inclusion of CHEMBL1163125 measurements of binding assays (type B) with IC50 values in nM
- Selection of exact (=) measurements only; < and > values were excluded as they provide no exact potency value; < values originated almost entirely from patent potency categories (e.g. < 100 nM)
- Exclusion of transcription errors and ChEMBL-flagged duplicates. 'Outside typical range' measurements were kept due to inconsistent application of the flag throughout the data
- Exclusion of (−)-JQ1 (CHEMBL3220211), the inactive enantiomer used as a negative control, due to contradictory values (likely mis-assigned stereochemistry)
- Selection of primary sources (scientific literature and patent data); secondary compilations (EUbOPEN, SGC, chemical probe sets) were excluded, as they partly re-report previously published values and contributed < 60 measurements
- Exclusion of patent values at 550, 2,750, 5,500 and 50,500 nM, identified as midpoints of IC50 ranges recorded as exact values
- Extraction of relevant columns to minimise excessive data, molecular weight refers to the parent compound without salt counterions

### Transformation: Python, `scripts/BRD4_IC50-conversion.py`
- Calculation of pIC50 from standard values, as pChEMBL values are missing for measurements with validity comments
- Reclassification of assay format (derived from ChEMBL's BAO label, corrected for protein expression hosts Sf9, BL21 and CHO) and annotation of assay technology (Alpha, TR-FRET, FP) 
- Separation of mutant BRD4 measurements from the regular BRD4 data as mutations may alter the binding site and potencies are not directly comparable to wild-type BRD4
- Condensation of protein assay-based measurements into aggregates per molecule: cell-, tissue- and organism-based values were excluded due to permeability and downstream effects on IC50 results; key information is the median pIC50, aggregation on the logarithmic scale weights fold-differences equally and the median limits the influence of outliers
- Classification into active and inactive inhibitors: threshold of pIC50 >= 6 (IC50 <= 1 µM), a common cutoff for active compounds

## Results
[![BRD4 dashboard](https://public.tableau.com/static/images/br/brd4_analysis/Dashboard/1.png)](https://public.tableau.com/views/brd4_analysis/Dashboard)

Interactive version: [Tableau Public](https://public.tableau.com/views/brd4_analysis/Dashboard)

Molecules without Rule-of-5 data in ChEMBL (n = 28) are not shown in the Rule-of-5 charts.

## Limitations
- Measurements on different BRD4 constructs (BD1, BD2, tandem) were combined, as the domain could not be reliably determined from assay descriptions; for domain-selective compounds, the median pIC50 therefore represents neither domain (the largest spreads in the data result from BD2-selective inhibitors)
- Differences in assay conditions between laboratories may be reflected in the spread of repeat measurements of the same molecule
- Actives outnumber inactives (3,879 vs. 1,990), likely reflecting publication bias towards optimised compound series
- Stereochemical annotation may be unreliable, so that measurements of one enantiomer may be assigned to the other

## How to run
1. Clone the repository and create the virtual environment:
    ```
    python -m venv venv
    .\venv\Scripts\activate 
    pip install -r requirements.txt
    ```

2. Download the ChEMBL SQLite database (`chembl_37_sqlite.tar.gz`) from the
   [ChEMBL FTP site](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/), unpack it (`tar -xzf chembl_37_sqlite.tar.gz`) and move `chembl_37.db` to `data/raw/`.

3. Run the data conversion script:
    ```
    python scripts/BRD4_IC50-conversion.py
    ```

The output is written to `data/processed/`:
- `brd4_measurements.csv`: all cleaned exact measurements and calculation of pIC50
- `brd4_mutants.csv`: separated data on mutant BRD4
- `brd4_molecules.csv`: one row per molecule including activity label
