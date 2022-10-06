#Rscript data.R

gcc -o3 rj-seed.c -lm -o rj-seed.ce -g

# as previous
./rj-seed.ce newest-data-old-mito-WILD.csv 1 50 100 0 100000 > tmpn1.1 &
./rj-seed.ce newest-data-old-mito-WILD.csv 2 50 100 0 100000 > tmpn1.2 &
./rj-seed.ce newest-data-old-mito-MSH1.csv 1 50 100 0 100000 > tmpn2.1 &
./rj-seed.ce newest-data-old-mito-MSH1.csv 2 50 100 0 100000 > tmpn2.2 &
./rj-seed.ce newest-data-old-plastid-MSH1.csv 1 50 100 0 100000 > tmpn3.1 &
./rj-seed.ce newest-data-old-plastid-MSH1.csv 2 50 100 0 100000 > tmpn3.2 &

# just new data, to compare old predictions
./rj-seed.ce newest-data-new-mito-WILD.csv 1 50 100 0 100000 > tmpn4.1 &
./rj-seed.ce newest-data-new-mito-WILD.csv 2 50 100 0 100000 > tmpn4.2 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 1 50 100 0 100000 > tmpn5.1 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 2 50 100 0 100000 > tmpn5.2 &
# just new data, enforcing linear dev model
./rj-seed.ce newest-data-new-mito-WILD.csv 1 50 100 -1 100000 > tmpn6.1 &
./rj-seed.ce newest-data-new-mito-WILD.csv 2 50 100 -1 100000 > tmpn6.2 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 1 50 100 -1 100000 > tmpn7.1 &
./rj-seed.ce newest-data-new-plastid-MSH1.csv 2 50 100 -1 100000 > tmpn7.2 &

# all data together
./rj-seed.ce newest-data-all-mito-WILD.csv 1 50 100 0 100000 > tmpn8.1 &
./rj-seed.ce newest-data-all-mito-WILD.csv 2 50 100 0 100000 > tmpn8.2 &
./rj-seed.ce newest-data-all-plastid-MSH1.csv 1 50 100 0 100000 > tmpn9.1 &
./rj-seed.ce newest-data-all-plastid-MSH1.csv 2 50 100 0 100000 > tmpn9.2 &

# test beds
#./rj-seed.ce data-bed-1.csv 50 100 0 100000 > tmp1c &
#./rj-seed.ce data-bed-2.csv 50 100 0 100000 > tmp2c &
#./rj-seed.ce data-bed-3.csv 50 100 0 100000 > tmp3c &
#./rj-seed.ce data-bed-4.csv 50 100 0 100000 > tmp4c &
#./rj-seed.ce data-bed-5.csv 50 100 0 100000 > tmp5c &
#./rj-seed.ce data-bed-6.csv 50 100 0 100000 > tmp6c &
#./rj-seed.ce data-bed-7.csv 50 100 0 100000 > tmp7c &


