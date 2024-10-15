# Peptide-arrays-for-Chagas-disease


## Linux (Bash ejecutable)
./script.R
Ej: ./03_calculate_peaks.R --main_folder /ruta/a/folder --testing FALSE --sources AR,BO,BR --min_amount_of_peptides_in_peak 3


## Linux (Bash)
Rscript 01_pools_normalize_data.R
Rscript 01_pools_normalize_data.R "C:/Users/Ramiro/Documents/GitHub/Peptide-arrays-for-Chagas-disease" TRUE "AR,BO,BR" 


## Para Windows (PowerShell)
& "C:\Program Files\R\R-4.2.1\bin\Rscript.exe" "01_pools_normalize_data.R"
& "C:\Program Files\R\R-4.2.1\bin\Rscript.exe" "01_pools_normalize_data_terminal.R" "C:/Users/Ramiro/Documents/GitHub/Peptide-arrays-for-Chagas-disease" TRUE "AR,BO,BR" 