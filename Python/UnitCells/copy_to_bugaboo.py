import numpy
import os
import subprocess

def copy_files():
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    units=[1,4,6,8,10]
    scale=[1,2,4,8,10,20]
    scale_str=[str(x) for x in scale]
    subprocess.call(['scp','jos_parallel.cpp','shurik@bugaboo.westgrid.ca:/home/shurik/UnitCells/'])
    subprocess.call(['scp','govern.py','shurik@bugaboo.westgrid.ca:/home/shurik/UnitCells/'])
    subprocess.call(['scp','binary.pbs','shurik@bugaboo.westgrid.ca:/home/shurik/UnitCells/'])
    for counter,dir_temp in enumerate(capillary_str):
        os.chdir(dir_temp)
        subprocess.call(['scp','geometry.dat','shurik@bugaboo.westgrid.ca:/home/shurik/UnitCells/'+dir_temp+"/"])
        subprocess.call(['scp','velx0200000.dat','shurik@bugaboo.westgrid.ca:/home/shurik/UnitCells/'+dir_temp+"/"])
        subprocess.call(['scp','vely0200000.dat','shurik@bugaboo.westgrid.ca:/home/shurik/UnitCells/'+dir_temp+"/"])
        print "Done with ",dir_temp 
        os.chdir("..")

if __name__=="__main__":
    #modify_file()
    #run_simulations()
    #copy_files("42",5,4)
    #copy_files("60",10,5)
    copy_files()
