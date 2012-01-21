import numpy
import os
import subprocess

def copy_files():
    capillary_str=[str(x) for x in [21,42,60,84]]
    units=[1,4,6,8,10]
    scale=[1,2,4,8,10,20]
    scale_str=[str(x) for x in scale]
    subprocess.call(['scp','jos_parallel.cpp','shurik@checkers.westgrid.ca:/home/shurik/UnitCellsOther/'])
    subprocess.call(['scp','govern.py','shurik@checkers.westgrid.ca:/home/shurik/UnitCellsOther/'])
    subprocess.call(['scp','binary.pbs','shurik@checkers.westgrid.ca:/home/shurik/UnitCellsOther/'])
    for counter,dir_temp in enumerate(capillary_str):
        os.chdir(dir_temp)
        subprocess.call(['scp','geometry.dat','shurik@checkers.westgrid.ca:/home/shurik/UnitCellsOther/'+dir_temp+"/"])
        subprocess.call(['scp','velx0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/UnitCellsOther/'+dir_temp+"/"])
        subprocess.call(['scp','vely0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/UnitCellsOther/'+dir_temp+"/"])
        print "Done with ",dir_temp 
        os.chdir("..")

if __name__=="__main__":
    #modify_file()
    #run_simulations()
    #copy_files("42",5,4)
    #copy_files("60",10,5)
    copy_files()
