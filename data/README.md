# BioMarker Pilot Dataset - README

## Purpose of the Dataset

This data set was created as the first pilot deliverable of the BioMarker Project (BioMarker: Development of an Integrated Software Platform for Automated Analysis, Structural Mapping, and Functional Prediction of Protein Sequences Using a Bioinformatics Approach, TÜBİTAK 2209-A)

The purpose of this pilot dataset is not to train a machine learning model. It was built to test and validate the computations and error-handling logic of the BioMarker software (sequence input, physicochemical analysis, hydropathy profiling, secondary structure prediction, motif scanning and BLAST-based homology analysis modules)

The pilot set contains 20 real, manually selected UniProt protein records, balanced across six structural/functional categories (soluble, membrane, secreted, nuclear, mitochondrial, and plant/fungal/bacterial proteins)

## Data Sources Used

All protein records were retrieved from official, publicly accessible bioinformatics databases:

|SOURCE|DATA RETRIEVED|
|-|-|
|UniProtKB/Swiss-Prot|Protein sequence, protein name, gene name, organism, taxonomy ID, sequence length, subcellular location, reviewed status|
|UniProt REST API (rest.uniprot.org)|Programmatic retrieval of the fields above in JSON format|
|RCSB PDB|3D structures, PDB IDs|
|PROSITE|Conserved motif and domain information contained of proteins, PROSITE ID|

Only Reviewed (Swiss-Prot) entries were used. Source URLs and retrieval dates for every record are logged in ‘source\_manifest.csv’

## Protein Selection and Exclusion Criteria

Full criteria are documented in ‘selection\_criteria.md’. Insummary:

**Included:**

* Reviewed (Swiss-Prot) entries only
* Complete amino acid sequences
* Valid, accessible UniProt accession numbers
* Clear protein name and organism information
* A balanced representation of protein categories and lengths

**Excluded (but logged, not deleted):**

* Fragment-flagged sequences
* Empty, truncated, or unverifiable sequences
* Duplicate accessions or duplicate sequences under different names
* Sequences with excessive ambiguous amino acid characters (e.g. X, Z)
* Records that could not be retrieved from the source database

Excluded records and their exclusion reasons are kept in ‘excluded\_records.csv’

## File Descriptions

|FILE|DESCRIPTION|
|-|-|
|data/raw/biomarker\_proteins\_20.xlsx and data/raw/biomarker\_proteins\_20.csv|Main data table containing all 20 proteins with full meta data (BioMarker\_ID, UniProt\_Accession, Protein\_Name, Organism, Sequence, QC\_Status, etc.)|
|data/sequences/biomarker\_sequences\_20.fasta|All 20 protein sequences in FASTA format. Header format: >BioMarker\_ID|
|docs/source\_research.xlsx|General research notes on each official data source (UniProt, PDB, PROSITE, InterPro, BLAST): data types, Access methods, licensing|
|data/manifest/source\_manifest.csv|Per-record audit trail — which database, URL, retrieval method, retrieval date, and WHO retrieved each record|
|docs/selection\_criteria.md|Written inclusion/exclusion rules applied when building the dataset|
|docs/data\_dictionary.xlsx|Column-by-column definitions for the main data table|
|data/excluded/excluded\_records.csv|Records considered but excluded, with the reason for exclusion|
|scripts/uniprot\_download.py|Python script that downloads protein data from UniProt and generates the CSV, FASTA, and manifest files|
|data/README.md|This file|

## How to Run the Code

a) Install the required package:
“pipinstallrequests”

b) Open ‘uniprot\_download.py’ and edit the ‘ACCESSIONS’ and ‘BIOMARKER\_IDS’ lists at the top of the file with the desired UniProt accession numbers and matching BioMarkerIDs (already filled in with the 20 selected proteins for this pilot set)

c) Run the script:
“python uniprot\_download.py”

d) The script will automatically generate:

* ‘biomarker\_proteins\_20.csv’ (main table)
* ‘biomarker\_sequences\_20.fasta’ (sequences)
* ‘source\_manifest.csv’ (retrieval log)
* ‘download\_errors.csv’ (only created if any download fails)

## Python and Package Versions Used

* Python: 3.10
* requests: 2.32.5
* charset\_normalizer: 3.4.4
* idna: 3.11
* urllib3: 2.6.3
* certifi: 2026.1.4

## Creation and Update Date

* **Created:** August 2026
* **Last updated:** August 21, 2026
* **Prepared by:** Şevvalnur ALTINOK, Tuba BİRSİN
* **Supervisor:** Dr. Öğr. Üyesi Tülay TURAN
* **Institution:** Burdur Mehmet Akif Ersoy University, Faculty of Engineering and Architecture, Department of Computer Engineering

## UPDATE HISTORY

August 11, 2026 - August 14, 2026: The first pilot dataset was created (20 proteins, base line columns filled)
August 11, 2026: Before starting protein selection, the criterion of preferring records with PDB and PROSITE information was determined.



## The Contents of The /data file



The pilot dataset (20 proteins) is stored here.

* 'raw/' — BioMarker\_Proteins\_20.csv / .xlsx (UniProt ID, gene name, family, PDB ID, PROSITE ID, sequence, weight, length, etc.)
* 'sequences/' — biomarker\_sequences\_20.fasta
* 'excluded/' — excluded\_records.csv (proteins excluded for not meeting selection criteria)
* 'manifest/' — source\_manifest.csv (TR/EN) — a record of the source from which each piece of data originates



## The Contents of The /docs file



This folder contains the methodology and documentation files for the Project.

* 'selection\_criteria.md' — Selection criteria for the 20 proteins
* 'data\_dictionary.xlsx' — Description of each field in the dataset
* 'source\_research.xlsx' — A file containing comprehensive information about the tools to be used in the project.
* 'cross\_check\_notes.md' — A document proving that the partners control the proteins they have selected.



## The Contents of The /scripts file

One-off or helper scripts are kept here (tools that are not part of the package but are used for data generation)

* 'uniprot\_download.py' — Script for fetching protein metadata and FASTA sequences from the UniProt REST API.



## The Contents of The /notebooks file



## The Contents of The /src file



## The Contents of The /tests file







