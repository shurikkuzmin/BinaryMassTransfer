import numpy
import os
import subprocess

def copy_files():
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    #subprocess.call(['scp','mass_periodic.cpp','shurik@checkers.westgrid.ca:/home/shurik/SamePeclet'])
    #subprocess.call(['scp','govern.py','shurik@checkers.westgrid.ca:/home/shurik/SamePeclet/'])
    #subprocess.call(['scp','binary.pbs','shurik@checkers.westgrid.ca:/home/shurik/SamePeclet/'])
    for counter,dir_temp in enumerate(capillary_str):
        os.chdir(dir_temp)
        subprocess.call(['scp','geometry.dat','shurik@checkers.westgrid.ca:/home/shurik/SamePeclet/'+dir_temp+"/"])
        subprocess.call(['scp','velx0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/SamePeclet/'+dir_temp+"/"])
        subprocess.call(['scp','vely0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/SamePeclet/'+dir_temp+"/"])
        print "Done with ",dir_temp 
        os.chdir("..")

if __name__=="__main__":
    #modify_file()
    #run_simulations()
    #copy_files("42",5,4)
    #copy_files("60",10,5)
    copy_files()
