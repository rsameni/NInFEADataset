# after generating the .dat and .hea files, run this line in the data folder to read back the data for comparison with the original data files
for i in {1..60}; do rdsamp -r $i -p > foo$i; done;