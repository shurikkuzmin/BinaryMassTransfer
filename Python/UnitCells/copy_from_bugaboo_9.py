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
        #subprocess.call(['scp','geometry.dat','shurik@checkers.westgrid.ca:/home/shurik/UnitCells/'+dir_temp+"/"])
        #subprocess.call(['scp','velx0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/UnitCells/'+dir_temp+"/"])
        #subprocess.call(['scp','vely0200000.dat','shurik@checkers.westgrid.ca:/home/shurik/UnitCells/'+dir_temp+"/"])
        print "Done with ",dir_temp 
        os.chdir("..")

def copy_one_directory(dir_name):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    units=[4,6,8,10]
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    scales=[0.3,0.5,1,2,4,6,8,10,15,20]    
    scale_str=["03","05","1","2","4","6","8","10","15","20"]
    os.chdir("Bugaboo/"+dir_name)
    flag=False
    counter=capillary_str.index(dir_name)
    for scale_counter,scale in enumerate(scales):
        if velocities[counter]*scale>0.12:
            scale3=scales[scale_counter-1]
            scale2=scales[scale_counter-2]
            scale1=scales[scale_counter-3]
            flag=True 
            break
    if not flag:
        scale1=scales[-3]
        scale2=scales[-2]
        scale3=scales[-1] 
    print scale1,scale2,scale3,velocities[counter] 
    for scale in [scale1,scale2,scale3]:
        dir_name_scale=scale_str[scales.index(scale)]
        subprocess.call(['mkdir','-p',dir_name_scale])
        os.chdir(dir_name_scale)
        for unit in units:
            subprocess.call(['mkdir','-p',str(unit)])
            os.chdir(str(unit))            
            subprocess.call(['scp','shurik@bugaboo.westgrid.ca:/home/shurik/UnitCells9/'+dir_name+"/"+dir_name_scale+"/"+str(unit)+"/density*.dat","."])
            print "Done with ",unit
            os.chdir("..")
        os.chdir("..")
        print "Scale ",scale," is finished"
 



if __name__=="__main__":
    #modify_file()
    #run_simulations()
    #copy_files("42",5,4)
    #copy_files("60",10,5)
    #copy_files()
    copy_one_directory("9")
