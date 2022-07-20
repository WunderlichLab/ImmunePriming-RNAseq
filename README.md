# Immune Priming RNA-seq

## General Description:
Downstream analysis for immune priming RNA-seq. This includes code to generate all relevant normalized data, graphs, and carry out analysis to identify priming-specific, immune loitering, and potentiated recall genes.

## File Descriptions:
- **ImmunePriming_RNA-seq_Downstream.R:** R Script used to generate all analyses after count generation. All package dependencies are listed at the beginning of the file. Analysis was run using R Version 4.2.0.

- **SurvivalCode_Basic.R:** R script used to generate all survival curves and test for significantly different survival between conditions. All package dependencies are listed at the beginning of the file. Analysis was run using R Version 4.2.0.

- Files for Automated Dilution Plating Analysis:
  - **saveanalyze_macros.ijm:** ImageJ macro used to count individual colonies on a JPG image of a bacterial plate. To use follow these steps:
    - Open up ImageJ
    - Go to  Plugins > Macros > Install
    - Select the file
    - The macro names and hotkeys will be listed under Plugin > Macros
  - **organizer_cfu.py**: Adds metadata and exports colony counts to Excel. Requires installation of Pandas for your Python environment. 
  - **DilutionPlating.R:** R script for creating bacterial load plots over time using the automated outputs from **organizer_cfu.py**
