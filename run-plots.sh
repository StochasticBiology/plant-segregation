# plotting scripts -- !!! these won't work until the inference processes have been launched (run-inference.sh) and won't take a final form until those runs have completed.

# plot summaries of the experimental data
Rscript plot-data.R

# plot summaries and posteriors for the test bed data
Rscript plot-posteriors-test.R

# plot posteriors for the experimental data
Rscript plot-posteriors-data.R

# plot figures for the gene conversion simulations
Rscript sim-seg.R
