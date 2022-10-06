Rscript parse-newest.R
Rscript plot-data-newest.R

gcc -o3 rj-seed.c -lm -o rj-seed.ce -g

mkdir tmp

# as previous
./rj-seed.ce newest-data-old-mito-WILD.csv 1 50 100 0 100000 > tmp/tmpn1.1 &
./rj-seed.ce newest-data-old-mito-WILD.csv 2 50 100 0 100000 > tmp/tmpn1.2 &
./rj-seed.ce newest-data-old-mito-MSH1.csv 1 50 100 0 100000 > tmp/tmpn2.1 &
./rj-seed.ce newest-data-old-mito-MSH1.csv 2 50 100 0 100000 > tmp/tmpn2.2 &
./rj-seed.ce newest-data-old-plastid-MSH1.csv 1 50 100 0 100000 > tmp/tmpn3.1 &
./rj-seed.ce newest-data-old-plastid-MSH1.csv 2 50 100 0 100000 > tmp/tmpn3.2 &

# just new data, to compare old predictions
./rj-seed.ce newest-data-new-mito-WILD.csv 1 50 100 0 100000 > tmp/tmpn4.1 &
./rj-seed.ce newest-data-new-mito-WILD.csv 2 50 100 0 100000 > tmp/tmpn4.2 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 1 50 100 0 100000 > tmp/tmpn5.1 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 2 50 100 0 100000 > tmp/tmpn5.2 &
# just new data, enforcing linear dev model
./rj-seed.ce newest-data-new-mito-WILD.csv 1 50 100 -1 100000 > tmp/tmpn6.1 &
./rj-seed.ce newest-data-new-mito-WILD.csv 2 50 100 -1 100000 > tmp/tmpn6.2 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 1 50 100 -1 100000 > tmp/tmpn7.1 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 2 50 100 -1 100000 > tmp/tmpn7.2 &

# all data together
./rj-seed.ce newest-data-all-mito-WILD.csv 1 50 100 0 100000 > tmp/tmpn8.1 &
./rj-seed.ce newest-data-all-mito-WILD.csv 2 50 100 0 100000 > tmp/tmpn8.2 &
./rj-seed.ce newest-data-all-plastid-MSH1.csv 1 50 100 0 100000 > tmp/tmpn9.1 &
./rj-seed.ce newest-data-all-plastid-MSH1.csv 2 50 100 0 100000 > tmp/tmpn9.2 &

# test beds
#./rj-seed.ce Data/data-bed-1.csv 1 50 100 0 100000 > tmp/tmpb.1 &
#./rj-seed.ce Data/data-bed-2.csv 1 50 100 0 100000 > tmp/tmpb.2 &
#./rj-seed.ce Data/data-bed-3.csv 1 50 100 0 100000 > tmp/tmpb.3 &
#./rj-seed.ce Data/data-bed-4.csv 1 50 100 0 100000 > tmp/tmpb.4 &
#./rj-seed.ce Data/data-bed-5.csv 1 50 100 0 100000 > tmp/tmpb.5 &
#./rj-seed.ce Data/data-bed-6.csv 1 50 100 0 100000 > tmp/tmpb.6 &
#./rj-seed.ce Data/data-bed-7.csv 1 50 100 0 100000 > tmp/tmpb.7 &


