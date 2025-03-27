# Code for Analysis of Peptide Array Data

Software and code accompanying an upcoming book chapter (citation below). In this chapter, we provide annotated code for analyzing peptide microarray data from the CHAGASTOPE project (aka the Chagas Antigen and Epitope Atlas). The same type of analysis can be applied to other peptide microarray data, with adjustments to the script codes if there is a different array design.

The repository contains a series of scripts for performing the various steps involved in the analysis and visualization of CHAGASTOPE data. The scripts are sequentially organized in code folder, with the order indicated by the number in the script name (e.g., 01_pools_normalize_data). The sequential steps include normalization of array data across samples, smoothing to remove outliers, calculation/detection of reactive (antigenic) peaks and regions, and analysis of single-residue mutagenesis scanning (Alanine Scans) of epitopes as described in (Ricci et al. 2023, DOI: [10.1038/s41467-023-37522-9](https://10.1038/s41467-023-37522-9)).

The code is divided into two sections: all main scripts contain the primary code, which reads and processes the parameters called and then run the main function. They are essentially wrappers for the main function. This main function, along with other auxiliary functions are  located in the functions sub-folder. Each script is numbered consistently in both sections. This separation is intended for clarity: the main script is concise, allowing users to quickly view the parameters in use and is designed for those who do not wish to modify the code. In contrast, the scripts containing all auxiliary functions are intended for users who wish to delve into the details of the process and make modifications. Parameters are located in the config or parameters sections of the respective scripts.

Both the main scripts and those in the functions folder must be downloaded for the code to function correctly.
The code can be executed either through the terminal (Bash or Windows PowerShell), which allows for modifying arguments passed to the script, or by running the script in a number of development environments that support the R programming language (RStudio is recommended).

A minimal set of data for testing purposes is provided in this repository in the `data` folder. The complete set of data from the Chagas Antigen and Epitope Atlas is available from the ArrayExpress database under accession numbers: E-MTAB-11651 and E-MTAB-11655.

## INSTALL 

All code has been tested with R version 4.3.3 and Bioconductor v3.18. 

### Instructions for Ubuntu Linux 22.04

The default R version in Ubuntu-24.04 is 4.3.3. This has been tested. The repository has [renv](https://rstudio.github.io/renv/articles/renv.html) config files that will auto install all required dependencies once an initial interactive R session is opened from the `code` directory. This interactive R session is only required once. 

```
# in the linux terminal

# update package indices
sudo apt update

# install R 
sudo apt install r-base-core libx11-dev

# clone this repo
git clone https://github.com/trypanosomatics/Peptide-arrays-for-Chagas-disease.git

# enter the code directory and make all R scripts executable
cd Peptide-arrays-for-Chagas-disease/code/
chmod +x *.R

# install all dependencies
# to start auto-installation, please initiate an interactive R session
R
```

### Instructions for Windows

Download and install [R-4.3.3](https://cran.r-project.org/bin/windows/base/old/4.3.3/), and maybe also [Rtools43](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) (which may be required to build R packages from source). 
Also, download and install [git for Windows](https://git-scm.com/downloads/win). Note that the R executables are not by default added to the PATH (see the [R for Windows FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Rcmd-is-not-found-in-my-PATH_0021)).

Then in a powershell terminal: 

```
# clone this repo
git clone https://github.com/trypanosomatics/Peptide-arrays-for-Chagas-disease.git

# enter the code directory 
cd Peptide-arrays-for-Chagas-disease/code/

# install all dependencies
# to start auto-installation, please initiate an interactive R session
R.exe 
```

### Install all dependencies (Linux or Windows)

For both Linux or Windows, the instructions are the same. Within an initial interactive R session, renv should bootstrap itself. You will notice a message saying "*One or more packages recorded in the lockfile are not installed.*" After renv is installed and ready, do: 

```
# this is all run within an interactive R session

# install all dependencies using renv::restore
# this will install all requirements listed in the renv.lock file
# and respecting all configs in renv/settings.json
# you should respond to occasional questions
renv::restore()

# check that everything is OK
renv::status()

# if you get this, you're OK to go!
No issues found -- the project is in a consistent state.

```

## USAGE

The more detailed descriptions of what each script does are written in the book chapter, but briefly, all code can be tested easily with a minimal data set provided in this repo. Some examples shown below: 

### Example Bash (Linux)

```
# minimal, assumes all testing data and defaults parameters
./01_pools_normalize_data.R --testing TRUE

# custom main folder and only some specific datasets
./01_pools_normalize_data.R --main_folder /path/to/folder --sources AR,BO,BR
```

### Example Windows PowerShell

For Windows / Poweshell scripts are passed as arguments to the `Rscript.exe` executable. The examples follow the same spirit as above. 

```
Rscript.exe '01_pools_normalize_data.R' --testing TRUE
Rscript.exe '01_pools_normalize_data.R' --main_folder "/path/to/folder" --testing FALSE --sources "AR,BO,BR"
```

## Citation 

If you use this code, please cite: 

Guadalupe Romer, Ramiro B Quinteros, Fernán Agüero. Software and tools to analyze high-density peptide array data for the Chagas Antigen and Epitope Atlas (2025). In: Trypanosoma cruzi infection: methods and protocols (Karina A Gómez & Carlos A Busgaglia, eds), Methods in Molecular Biology (series). Springer / Humana Press. **In process.**

Bibtex citation (will be updated soon): 
```
@incollection{romer_25_software,
  author    = {Romer G, Quinteros RB, Agüero F},
  title     = {Software and tools to analyze high-density peptide array data for the Chagas Antigen and Epitope Atlas},
  year      = {2025},
  chapter   = {},
  pages     = {},
  editor    = {Gómez KA, Buscaglia CA},
  booktitle = {Trypanosoma cruzi infection: methods and protocols},
  series    = {Methods in Molecular Biology},
  publisher = {Humana Press},
  volume    = {},
  number    = {},
  issn      = {},
  doi       = {},
  url       = {},
}
```
