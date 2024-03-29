# plant-segregation
Code to infer dynamics of oDNA segregation in plants from experimental data

Preprint here https://www.biorxiv.org/content/10.1101/2022.11.07.515340v2

Requirements
------
R with packages `ggplot2`, `gridExtra`, `readxl`, and `stringr`; C compiler. The wrapper script `run-inference.sh` is for a Bash environment: if you are running in a different environment, take a look at the commands there to see what needs doing.

Outline
------
The pipeline consists of R code to parse experimental datafiles and synthetic test cases (in `Data/`), C code `rj-inference.c` to perform RJMCMC inference on these parsed outputs, C code `sim-seg.c` to simulate gene conversion, and R code for plotting all outputs. The wrapper script `run-inference.sh` runs the parser then launches all the simulations. The RJMCMC code takes a while -- there are 5 test cases and 14 data cases for analysis, and the larger ones take perhaps a few days on a modern processor. The script currently forks each process so they are run in parallel.

The target of the inference processes is the various extents of segregation (increasing heteroplasmy variance) through development and between generations, and (in parallel) the structure of the developmental model linking these observations.

For computational efficiency, the RJMCMC code precomputes a lookup table of probability density values for the core "Kimura" distribution involved in the inference process. This precomputation step takes some minutes, but when it is completed a file `dbtrj-50-100.csv` is placed in the working directory where it can be quickly read by other processes. A downside to running the code via the wrapper script is that all 19 cases will start at the same time, none will see an existing version of this file, and all will spend time precomputing it. To get around this, we include a zipped version of this precomputed file `dbtrj-50-100.csv.gz` in the repo; the wrapper script unzips it before running code. If the parameters of the simulation are changed, a new precomputation will be required -- to avoid the inefficiency above, one can run a single case, wait until the file appears, then set the other instances running.

Each inference run produces an output CSV files storing posterior samples. These are used by followup plot scripts to produce figures. The wrapper `run-plots.R` calls each plotting script in turn to produce these. 

Update
------
As of March 2023 the repo has been updated with code `plot-posteriors-data-postrev.R` and `run-pruned.sh` which remove an outlying lineage from the dataset and re-run the analysis, producing more plots to describe the system without this lineage. In July 2023 various other `-postrev.R` plotting scripts were added to make styling changes to the plots from peer review suggestions.
