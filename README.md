# plant-segregation

Code to infer dynamics of oDNA segregation in plants from experimental data.

Requirements
------
R with packages `ggplot2`, `gridExtra`, `readxl`, and `stringr`; C compiler. The wrapper script `run-inference.sh` is for a Bash environment: if you are running in a different environment, take a look at the commands there to see what needs doing.

Outline
------
The pipeline consists of R code to parse experimental datafiles, C code `rj-inference.c` to perform RJMCMC inference on these parsed outputs, C code `sim-seg.c` to simulate gene conversion, and R code for plotting all outputs. The wrapper script `run-inference.sh` runs the parser then launches all the simulations. The RJMCMC code takes a while -- there are 5 test cases and 14 data cases for analysis, and the larger ones take perhaps a few days on a modern processor. The script currently forks each process so they are run in parallel.

Each inference run produces an output CSV files storing posterior samples. These are used by followup plot scripts to produce figures. The wrapper `run-plots.R` calls each plotting script in turn to produce these. 
