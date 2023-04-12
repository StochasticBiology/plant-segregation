# ensure the prelims from run-inference.sh have been run first

# remove particular highly-segregated family line from dataset
grep -v "104,4,[23]" newest-data-old-mito-MSH1.csv > newest-data-oldpruned-mito-MSH1.csv

# do inference on this pruned dataset
./rj-inference.ce newest-data-oldpruned-mito-MSH1.csv 1 50 100 0 100000 > tmp/tmpn10.1 &
./rj-inference.ce newest-data-oldpruned-mito-MSH1.csv 2 50 100 0 100000 > tmp/tmpn10.2 &
