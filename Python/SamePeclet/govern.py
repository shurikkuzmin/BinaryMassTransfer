#import numpy
import os
import subprocess
def modify_file(peclet_str):
    f=open("binary.pbs","r+")
    f2=open("binary_new.pbs","w")
    for counter,line in enumerate(f):
        if line.find("TOCHANGE")!=-1:
            line=line.replace("TOCHANGE",peclet_str)

        f2.write(line)
            
    f.close()
    f2.close()


def run_simulations():
    capillary_str=[str(x) for x in [9,21,42,60,84]]

    for capillary in capillary_str:
        dir_capillary_name=str(capillary)
        subprocess.call(['cp','main.out',dir_capillary_name+"/"])
        subprocess.call(['cp','binary.pbs',dir_capillary_name+"/"])
        os.chdir(dir_capillary_name)
        for peclet in range(2000,18000,2000):
            dir_name=str(peclet)
            subprocess.call(['mkdir','-p',dir_name])
            subprocess.call(['cp','main.out',dir_name+"/"])
            subprocess.call(['cp','binary.pbs',dir_name+"/"])
            subprocess.call(['cp','geometry.dat',dir_name+"/"])
            subprocess.call(['cp','velx0200000.dat',dir_name+"/"])
            subprocess.call(['cp','vely0200000.dat',dir_name+"/"])
            os.chdir(dir_name)
            modify_file(str(peclet))
            #subprocess.call(['qsub','binary_new.pbs'])
 
            os.chdir("..")
        os.chdir("..")

if __name__=="__main__":
    #modify_file()
    run_simulations()
