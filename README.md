# BioMarker

An integrated bioinformatics software platform for the automated analysis, structural mapping, and functional prediction of protein sequences.

It is being developed under the TÜBİTAK 2209-A Support Program for University Students' Research Projects in Türkiye.

## Project Summary

## Project Summary

Protein analyses today typically require a fragmented workflow: a web tool for physicochemical properties, another site for hydropathy, and different servers for motif scanning and homology analysis, etc.

BioMarker aims to consolidate all these analyses into a single open-source platform.

## Project Status (Last Update: August 2026)

* \[x] **Phase 0 — Pilot Veri Seti:** Metadata for 20 proteins, FASTA, selection criteria, and the source manifest have been completed.
* \[x] **Phase 1 (Week 1-2):** Literature review and technological infrastructure research
* \[x] **Phase 2 (Week 3-4):** Physicochemical profile + Hydropathy modules
* \[ ] **Phase 3 (Week 5-6):** Chou-Fasman secondary structure estimation
* \[ ] **Phase 4 (Week 7-9):** PyQt5 graphical user interface (GUI)
* \[ ] **Phase 5 (Week 10-12):** Motif scanning (PROSITE) + BLAST homology analysis
* \[ ] **Phase 6 (Week 13-15):** Interactive visualization (Plotly) + PDF reporting
* \[ ] **Phase 7 (Week 16-18):** 3D structure visualization (NGLviewer) integration
* \[ ] **Phase 8 (Week 19-24):** Performance testing, documentation, packaging

## Repo Structre


## Repo Structure

```text
biomarker/
├── data/ # Pilot dataset (protein metadata, FASTA, etc.)
├── docs/ # Selection criteria, data dictionary, methodology
├── src/biomarker/ # Main software package
│   ├── data_io/ # Data entry, validation, UniProt retrieval
│   ├── analysis/ # Physicochemical, hydropathy, secondary structure, motif, BLAST
│   ├── visualization/ # Static + interactive graphics, 3D structure visualization
│   ├── reporting/ # PDF report generation
│   └── gui/ # PyQt5 Interface
├── tests/ # Unit tests
└── notebooks/ # Exploratory data analysis
```


## Technologies Used

* **Language:** Python 3.x
* **Bioinformatics:** Biopython (SeqIO, SeqUtils, Entrez, BLAST)
* **Data Processing:** NumPy, Pandas
* **Visualization:** Matplotlib, Seaborn, Plotly
* **Interface:** PyQt5
* **Reporting:** WeasyPrint

## Setup

'''

bash
git clone https://github.com/sevvalAltinok16/biomarker.git
cd biomarker
pip install -r requirements.txt


'''

## About Dataset

The pilot dataset contains 20 proteins; for each, the UniProt ID, gene name, family, cellular localization, PDB ID, PROSITE ID, sequence, molecular weight, and more are recorded. 

Please refer to the `docs/data_dictionary.xlsx` and `docs/selection_criteria.md` files for details.

## Licence

This project is being developed as open source. See at 'LICENSE'

## Source

This project is being conducted within the scope of the TÜBİTAK 2209-A program at the Department of Computer Engineering, Faculty of Engineering and Architecture, Mehmet Akif Ersoy University.

