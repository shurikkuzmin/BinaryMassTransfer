#!/bin/bash

#mkdir -p 9 21 42 60 84
scp mass_periodic.cpp "shurik@checkers.westgrid.ca":/home/shurik/ScalingPeclet/$i/
scp govern.py "shurik@checkers.westgrid.ca":/home/shurik/ScalingPeclet/$i/
for i in 9 21 42 60 84
do
  cd $i
  scp geometry.dat "shurik@checkers.westgrid.ca":/home/shurik/ScalingPeclet/$i/ 
  scp vel[x,y]0200000.dat "shurik@checkers.westgrid.ca":/home/shurik/ScalingPeclet/$i/
  cd ..
  echo "Done with $i"
done
