import numpy
import os
import subprocess
def modify_file(scale_str):
    f=open("binary.pbs","r+")
    f2=open("binary_new.pbs","w")
    for counter,line in enumerate(f):
        if line.find("TOCHANGE")!=-1:
            line=line.replace("TOCHANGE",scale_str)
        f2.write(line)
            
    f.close()
    f2.close()


def run_simulations():
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    scale_str=[str(x) for x in [2,4,6,8,10,15,20,40]]
    for counter,dir_temp in enumerate(capillary_str):
        subprocess.call(['cp','main.out',dir_temp+"/"])
        subprocess.call(['cp','binary.pbs',dir_temp+"/"])
        os.chdir(dir_temp)
        for scale in range(0,len(scale_str)-counter):
            subprocess.call(['mkdir','-p',scale_str[scale]])
            subprocess.call(['cp','main.out',scale_str[scale]+"/"])
            subprocess.call(['cp','binary.pbs',scale_str[scale]+"/"])
           
            subprocess.call(['cp','geometry.dat',scale_str[scale]+"/"])
            subprocess.call(['cp','velx0200000.dat',scale_str[scale]+"/"])
            subprocess.call(['cp','vely0200000.dat',scale_str[scale]+"/"])
            os.chdir(scale_str[scale])
            modify_file(scale_str[scale])
            #subprocess.call(['qsub'],'binary_new.pbs']) 
            os.chdir("..")
            #subprocess.call(['qsub','binary_new.pbs'])
        os.chdir("..")

if __name__=="__main__":
    #modify_file()
    run_simulations()
