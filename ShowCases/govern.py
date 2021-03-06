#import numpy
import os
import subprocess
def modify_file(length,width_value,force):
    f=open("binary.pbs","r+")
    f2=open("binary_new.pbs","w")
    for counter,line in enumerate(f):
        if line.find("TOCHANGE")!=-1:
            line=line.replace("TOCHANGE1",str(length))
            line=line.replace("TOCHANGE2",str(width_value))
            line=line.replace("TOCHANGE3",str(force))
        f2.write(line)
            
            
    f.close()
    f2.close()
def run_simulations():
    capillary=[0.01*x for x in range(9,93,10)]
    capillary_str=[str(x) for x in range(9,93,10)]
    force_init=6e-6/16;
    width=0.1
    for length in range(200,1400,150):
        subprocess.call(['mkdir','-p',str(length)])
        subprocess.call(['cp','main.out',str(length)+"/"])
        subprocess.call(['cp','binary.pbs',str(length)+"/"])
        os.chdir(str(length))
        for i in range(0, len(capillary)):
            dir_temp=capillary_str[i]
            subprocess.call(['mkdir','-p',dir_temp])
            subprocess.call(['cp','main.out',dir_temp+"/"])
            subprocess.call(['cp','binary.pbs',dir_temp+"/"])
            os.chdir(dir_temp)
            subprocess.call(["mkdir","-p","tmp"])
            ratio=capillary[i]/0.05
            ratio2=1.0/(1.0+(length-200.0)*(length-200.0)/(150.0*150.0*26.0))
            force=force_init*ratio*ratio2
            width_value=int(width*200)
            modify_file(length,width_value,force)
            subprocess.call(['qsub','binary_new.pbs'])
        
            os.chdir("..")
        os.chdir("..")

if __name__=="__main__":
    #modify_file()
    run_simulations()
