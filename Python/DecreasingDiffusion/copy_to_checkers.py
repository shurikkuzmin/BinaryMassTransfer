import numpy
import os
import subprocess

def copy_files(dir_name):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    #units=[4,6,8]
    #scales=[0.3,0.5,1,2,4,6,8,10,15,20]    
    #scale_str=["03","05","1","2","4","6","8","10","15","20"]
    #os.chdir(dir_name)
    #scale_str=[str(x) for x in scale]
    subprocess.call(['scp','jos_parallel.cpp','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion/'])
    subprocess.call(['scp','mass_periodic.cpp','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion'])
    subprocess.call(['scp','govern.py','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion/'])
    subprocess.call(['scp','binary.pbs','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion/'])
    subprocess.call(['scp','binary_periodic.pbs','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion/'])
    os.chdir(dir_name)
    #for counter,dir_temp in enumerate(capillary_str):
    #    os.chdir(dir_temp)
    subprocess.call(['scp','geometry.dat','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion/'+dir_name+"/"])
    subprocess.call(['scp','velx0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion/'+dir_name+"/"])
    subprocess.call(['scp','vely0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/DecreasingDiffusion/'+dir_name+"/"])
    print "Done with ",dir_name 
    os.chdir("..")

if __name__=="__main__":
    #modify_file()
    #run_simulations()
    #copy_files("42",5,4)
    #copy_files("60",10,5)
    copy_files("9")
