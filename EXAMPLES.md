# Example commands and code 

These are all example commands and code as provided in the book chapter.
Because the book chapter has these code boxes rendered as images, we are
providing the examples as text here, to help readers copy and paste the
examples. 

The order and numbering of all examples follow the chapter sections. 

## 3. Methods

```bash
# this is a comment
command --parameter1 "argument 1" –-parameter2 "value1, value2"
```

### 3.1 Download the code from the GitHub Repository

```bash
git clone https://github.com/trypanosomatics/Peptide-arrays-for-Chagas-disease.git
```

### 3.2 Download the Data Set

Example code to download assay data using ‘wget’:

```bash
# change into the desired directory
cd Peptide-arrays-for-Chagas-disease
cd data/chagastope_data/inputs/02_pools_raw_data 
# and download the desired files
# only some examples are shown for brevity
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-11651/AR_NE_raw.tsv
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-11651/AR_PO_raw.tsv
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-11655/BO_E1_PO_raw.tsv
```

Example code to download the array mapping data using wget:

```bash
# change into the desired directory, for example for CHAGASTOPE-v1
cd data/chagastope_data/inputs/01_pools_array_design
# and get the data
wget https://figshare.com/ndownloader/files/38867358
# uncompress
unzip 38867358 
# Make sure the file name is: "Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv"
```

### 3.3 Make scripts executable

```bash
# from the top-level directory of the cloned repository
# change into the 'code' directory
cd Peptide-arrays-for-Chagas-disease/code
# and give the execute permission to all R script
chmod +x *.R
```

### 3.4 Normalize Peptide Array Pool Data

Default parameters: `./01_pools_normalize_data.R`


Custom parameters (example):
```bash
./01_pools_normalize_data.R \
  --main_folder /path/to/folder \
  --testing TRUE \
  --sources AR,BO,BR
```

### 3.5 Smooth Peptide Array Pool Data

Default parameters: `./02_pools_smooth_data.R`

Custom parameters (example):
```bash
./02_pools_smooth_data.R \
  --sources AR,BO,BR \
  --smoothing_mean_window_size 5 \
  --smoothing_median_window_size 1 \
  --smooth_borders_option "zeros" \
  --testing TRUE \
  --main_folder "/path/to/folder/" \
  --output_suffix "mean_1_median_5"
```

### 3.6 Analyze Antigenic Peaks from Peptide Array Data

Default parameters: `./03_calculate_peaks.R`

Custom parameters (example):
```bash
./03_calculate_peaks.R \
  --profile_data_suffix "mean_1_median_5" \
  --testing TRUE \
  --sources "AR,BO,BR" \
  --min_num_of_peptides_in_peak 2 \
  --sd_multiplier_for_cutoff 5
```

### 3.7 Analyze Antigenic Regions from Peptide Array Data

Default parameters: `./04_calculate_regions.R`

Custom parameters (example):
```bash
./04_calculate_regions.R \
  --peaks_tag "5SD_2pep" \
  --testing TRUE
```

### 3.8 Normalize Peptide Array Individual Sera Data

Default parameters: `./05_individual_serums_normalize_data.R`

Custom parameters (example):
```bash
./05_individual_serums_normalize_data.R \
  --main_folder /path/to/folder \
  --testing TRUE \
  --sources AR_P1,AR_E2,BR_P1,BR_E2,BO_P1,BO_E2
```

### 3.9 Smooth Peptide Array Individual Sera Data

Default parameters: `./06_individual_serums_smooth_data.R`

Custom parameters (example):
```bash
./06_individual_serums_smooth_data.R \
  --sources AR_P1,AR_E2,BR_P1 \
  --smoothing_median_window_size 5 \
  --smoothing_mean_window_size 1 \
  --smooth_borders_option "zeros" \
  --testing TRUE \
  --main_folder "/path/to/folder/" \
  --output_suffix "mean_1_median_5"
```

### 3.10 Analysis of single-residue mutagenesis data (alanine scans)

Default parameters: `./07_alanine_scan_analysis.R`

Custom parameters (example):
```bash
./07_alanine_scan_analysis.R \
  --main_folder /path/to/folder \
  --selected_protein "TCSYLVIO_004530" \
  --testing FALSE \
  --heatmap_plot TRUE \
  --sequence_logo_plot TRUE \
  --sequence_logo_source "all"
```

### 3.11 Plotting Smoothed Protein Profiles

#### For pooled samples: 

Default parameters: `./08_pools_plot_protein_profiles.R`

Custom parameters (example):
```bash
./08_pools_plot_protein_profiles.R \
  --sources "AR,BR" \
  --main_folder "/path/to/folder" \
  --profile_data_suffix "mean_5_median_3" \
  --sd_multiplier_for_cutoff 3 \
  --only_proteins_above 50 \
  --output_suffix "mean_5_median_3"
```

#### For individual samples:

Default parameters: `./09_individual_serums_plot_protein_profiles.R`

Custom parameters (example):
```bash
./09_individual_serums_plot_protein_profiles.R \
  --sources AR_P1, AR_P2, AR_P3, BO_P1, BO_P2, BR_P1, BR_P2 \
  --main_folder "/path/to/folder" \
  --protein "TCSYLVIO_004530" \
  --profile_data_suffix "mean_5_median_3" \
  --sd_multiplier_for_cutoff 3 \
  --output_suffix "mean_5_median_3"
```

#### For signal variations: 

Default parameters: `./10_plot_signal_variations.R`

Custom parameters (example):
```bash
./10_plot_signal_variations.R \
  --testing FALSE \
  --main_folder "/path/to/folder" \
  --type_data "individual_serums" \
  --smoothed_data_suffix "mean_5_median_3_smoothed_signals.tsv"
```

## Annex - Running the code in Windows

Go to the directory where you have the repository's code:

```powershell
cd "C:/Users/%USERPROFILE%/Documents/GitHub/Peptide-arrays-for-Chagas-disease/code"
```

### A1. Bootstrap the R installation and required packages

Run an interactive R session: `R.exe`

```R
# first let R bootstrap the installation of renv
# then run this
renv::restore()
# and when this is done check status with 
renv::status()
```

If you get: 'No issues found -- the project is in a consistent state.', you are
OK to go.

### 3.4 Normalize Peptide Array Pool Data

Default parameters (see NOTE 11): `Rscript '01_pools_normalize_data.R'`

Custom parameters:
```powershell
Rscript '01_pools_normalize_data.R' --main_folder "/path/to/folder" --testing FALSE --sources "AR,BO,BR"
```

### 3.5 Smooth Peptide Array Pool Data

Default parameters: `Rscript '02_pools_smooth_data.R'`

Custom parameters:
```powershell
 Rscript '02_pools_smooth_data.R' --sources AR,BO,BR --smoothing_median_window_size 5 --smoothing_mean_window_size 1 --smooth_borders_option "zeros" --testing FALSE --main_folder "/path/to/folder/" --output_suffix "mean_1_median_5"
```

### 3.6 Analyze Antigenic Peaks from Peptide Array Data

Default parameters: `Rscript '03_calculate_peaks.R'`

Custom parameters:
```powershell
Rscript '03_calculate_peaks.R' --sources AR,BO,BR --testing FALSE --main_folder "/path/to/folder/"
```

### 3.7 Analyze Antigenic Regions from Peptide Array Data

Default parameters: `Rscript '04_calculate_regions.R'`

Custom parameters:
```powershell
Rscript '04_calculate_regions.R' --testing FALSE --main_folder "/path/to/folder/"
```

### 3.8 Normalize Peptide Array Individual Sera Data

Default parameters: `Rscript '05_individual_serums_normalize_data.R'`

Custom parameters:
```powershell
Rscript '05_individual_serums_normalize_data.R' --sources AR_P1,AR_P2,AR_E1,AR_E2,BR_P1,BR_P2,BR_E1,BR_E2,BO_P1,BO_P2,BO_E1,BO_E2 --testing FALSE --main_folder "/path/to/folder/"
```

### 3.9 Smooth Peptide Array Individual Sera Data

Default parameters: `Rscript '06_individual_serums_smooth_data.R'`

Custom parameters:
```powershell
Rscript '06_individual_serums_smooth_data.R' --sources AR_P1,AR_P2,AR_E1,AR_E2,BR_P1,BR_P2,BR_E1,BR_E2,BO_P1,BO_P2,BO_E1,BO_E2 --testing FALSE --main_folder "/path/to/folder/" smoothing_median_window_size 4 --smoothing_mean_window_size 2 --smooth_borders_option "zeros"
```

### 3.10 Analysis of single-residue mutagenesis data (alanine scans)

Default parameters: `Rscript '07_alanine_scan_analysis.R'`

Custom parameters:
```powershell
Rscript '07_alanine_scan_analysis.R' --testing FALSE --main_folder "/path/to/folder/" --selected_protein TcCLB.511671.60
```

3.11 Plotting Smoothed Protein Profiles

Default parameters: `Rscript '08_pools_plot_protein_profiles.R'`

Custom parameters:
```powershell
Rscript '08_pools_plot_protein_profiles.R' --testing FALSE --main_folder "/path/to/folder/" --sources AR,BO,BR --protein  "TcCLB.511671.60"
```


