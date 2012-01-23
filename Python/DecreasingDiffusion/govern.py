import numpy
import os
import subprocess
def modify_file(scale_str,unit_str,scale_diff_str):
    f=open("binary.pbs","r+")
    f2=open("binary_new.pbs","w")
    for counter,line in enumerate(f):
        if line.find("TOCHANGE")!=-1:
            line=line.replace("TOCHANGE1",scale_str)
            line=line.replace("TOCHANGE2",unit_str)
            line=line.replace("TOCHANGE3",unit_str)
            line=line.replace("TOCHANGE4",scale_diff_str)
        f2.write(line)
            
    f.close()
    f2.close()

def modify_periodic_file(scale_str,unit_str,scale_diff_str):
    f=open("binary_periodic.pbs","r+")
    f2=open("binary_periodic_new.pbs","w")
    for counter,line in enumerate(f):
        if line.find("TOCHANGE")!=-1:
            line=line.replace("TOCHANGE1",scale_str)
            line=line.replace("TOCHANGE2",unit_str)
            line=line.replace("TOCHANGE3",unit_str)
            line=line.replace("TOCHANGE4",scale_diff_str)
        f2.write(line)
            
    f.close()
    f2.close()

def run_simulations(dir_name):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    units=[4,6,8,10]
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    scales=[0.3,0.5,1,2,4,6,8,10,15,20]    
    scale_str=["03","05","1","2","4","6","8","10","15","20"]
    scales_diff=[0.05,0.1,0.2,0.4,0.5,1.0]
    
    #for counter,dir_temp in enumerate(capillary_str):
    subprocess.call(['cp','main.out',dir_name+"/"])
    subprocess.call(['cp','main_diffusion.out',dir_name+"/"])
    subprocess.call(['cp','binary.pbs',dir_name+"/"])
    os.chdir(dir_name)
    
    flag=False
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
        dir_name=scale_str[scales.index(scale)]
        subprocess.call(['mkdir','-p',dir_name])
        subprocess.call(['cp','main.out',dir_name+"/"])
        subprocess.call(['cp','main_diffusion.out',dir_name+"/"])
        subprocess.call(['cp','binary.pbs',dir_name+"/"])
        subprocess.call(['cp','geometry.dat',dir_name+"/"])
        subprocess.call(['cp','velx0200000.dat',dir_name+"/"])
        subprocess.call(['cp','vely0200000.dat',dir_name+"/"])
        os.chdir(dir_name)
        for unit in units:
            subprocess.call(['mkdir','-p',str(unit)])
            subprocess.call(['cp','main.out',str(unit)+"/"])
            subprocess.call(['cp','main_diffusion.out',str(unit)+"/"])
            subprocess.call(['cp','binary.pbs',str(unit)+"/"])
            subprocess.call(['cp','geometry.dat',str(unit)+"/"])
            subprocess.call(['cp','velx0200000.dat',str(unit)+"/"])
            subprocess.call(['cp','vely0200000.dat',str(unit)+"/"])
            os.chdir(str(unit))            
            modify_file(str(scale),str(unit))
            subprocess.call(['qsub','binary_new.pbs']) 
            os.chdir("..")
            os.chdir("..")
            #subprocess.call(['qsub','binary_new.pbs'])
    os.chdir("..")

if __name__=="__main__":
    #modify_file()
    run_simulations("9")
