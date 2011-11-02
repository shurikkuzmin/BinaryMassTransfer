#!/bin/bash


for i in 3 5 8 10
do
  cd $i
  for j in `seq 1160000 280000 3960000`
  do
    scp "shurik@checkers.westgrid.ca":/home/shurik/MassSailfishJosDouble/$i/density_jos_par$j.dat .
  done
cd ..
  echo "Done with $i"
done
