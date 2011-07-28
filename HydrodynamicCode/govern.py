import numpy
import os
import subprocess
def modify_file(width_value,force):
    f=open("binary.pbs","r+")
    f2=open("binary_new.pbs","w")
    for counter,line in enumerate(f):
        if line.find("TOCHANGE")!=-1:
            line=line.replace("TOCHANGE1",str(width_value))
            line=line.replace("TOCHANGE2",str(force))
        f2.write(line)
            
            
    f.close()
    f2.close()
def run_simulations():
    capillary=numpy.linspace(0.03,0.9,30)
    capillary_str=[str(x) for x in range(3,93,3)]
    force_init=6e-6/16;
    width=numpy.linspace(0.04,0.17,30)
    for i in range(0, len(capillary)):
        dir_temp=capillary_str[i]
        os.mkdir(dir_temp)
        subprocess.call(['cp','main.out',dir_temp+"/"])
        subprocess.call(['cp','binary.pbs',dir_temp+"/"])
        os.chdir(dir_temp)
        os.mkdir("tmp")
        ratio=capillary[i]/0.05
        force=force_init*ratio
        width_value=int(width[i]*200)
        modify_file(width_value,force)
        #subprocess.call(['qsub','binary_new.pbs'])
        
        os.chdir("..")

if __name__=="__main__":
    #modify_file()
    run_simulations()
