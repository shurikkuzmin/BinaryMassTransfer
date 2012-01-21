#!/bin/bash

mkdir -p `seq 3 3 90`

for i in $(seq 3 3 90)
do
  cd $i
  scp "shurik@hermes.westgrid.ca":/home/shurik/HydroAlexandra/$i/tmp/velocityx300000.dat .
  scp "shurik@hermes.westgrid.ca":/home/shurik/HydroAlexandra/$i/tmp/velocityy300000.dat .
  scp "shurik@hermes.westgrid.ca":/home/shurik/HydroAlexandra/$i/tmp/phase300000.dat .
  cd ..

  echo "Done with $i"
done
