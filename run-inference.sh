# parse Excel data into CSV files of convenient for
Rscript parse-data.R

# plot summaries of these data
Rscript plot-data.R

# compile inference code
gcc -o3 rj-inference.c -lm -o rj-inference.ce -g

# make directory for incidental outputs (just for debugging and progress tracking really)
mkdir tmp/

# pull test bed data into this folder
cp Data/data-bed* .

# run inference on test beds
./rj-inference.ce data-bed-1.csv 1 50 100 0 100000 > tmp/tmpb.1 &
./rj-inference.ce data-bed-2.csv 1 50 100 0 100000 > tmp/tmpb.2 &
./rj-inference.ce data-bed-3.csv 1 50 100 0 100000 > tmp/tmpb.3 &
./rj-inference.ce data-bed-4.csv 1 50 100 0 100000 > tmp/tmpb.4 &
./rj-inference.ce data-bed-5.csv 1 50 100 0 100000 > tmp/tmpb.5 &

# run inference on existing datasets
./rj-inference.ce newest-data-old-mito-WILD.csv 1 50 100 0 100000 > tmp/tmpn1.1 &
./rj-inference.ce newest-data-old-mito-WILD.csv 2 50 100 0 100000 > tmp/tmpn1.2 &
./rj-inference.ce newest-data-old-mito-MSH1.csv 1 50 100 0 100000 > tmp/tmpn2.1 &
./rj-inference.ce newest-data-old-mito-MSH1.csv 2 50 100 0 100000 > tmp/tmpn2.2 &
./rj-inference.ce newest-data-old-plastid-MSH1.csv 1 50 100 0 100000 > tmp/tmpn3.1 &
./rj-inference.ce newest-data-old-plastid-MSH1.csv 2 50 100 0 100000 > tmp/tmpn3.2 &

# run inference on just new data, to compare old predictions
./rj-inference.ce newest-data-new-mito-WILD.csv 1 50 100 0 100000 > tmp/tmpn4.1 &
./rj-inference.ce newest-data-new-mito-WILD.csv 2 50 100 0 100000 > tmp/tmpn4.2 &
./rj-inference.ce newest-data-new-plastid-MSH1.csv 1 50 100 0 100000 > tmp/tmpn5.1 &
./rj-inference.ce newest-data-new-plastid-MSH1.csv 2 50 100 0 100000 > tmp/tmpn5.2 &

# run inference on all data together
./rj-inference.ce newest-data-all-mito-WILD.csv 1 50 100 0 100000 > tmp/tmpn8.1 &
./rj-inference.ce newest-data-all-mito-WILD.csv 2 50 100 0 100000 > tmp/tmpn8.2 &
./rj-inference.ce newest-data-all-plastid-MSH1.csv 1 50 100 0 100000 > tmp/tmpn9.1 &
./rj-inference.ce newest-data-all-plastid-MSH1.csv 2 50 100 0 100000 > tmp/tmpn9.2 &

