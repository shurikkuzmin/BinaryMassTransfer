#import numpy
import os
import subprocess

def run_simulations(main_dir_name):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    units=[4,6,8]
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    scales=[0.3,0.5,1,2,4,6,8,10,15,20]    
    scale_str=["03","05","1","2","4","6","8","10","15","20"]
    scales_diff=[2,5,10,20,40]
    counter=capillary_str.index(main_dir_name)
    
    #for counter,dir_temp in enumerate(capillary_str):
    os.chdir(main_dir_name)
    
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
    for scale in [scale2,scale3]:
        dir_name=scale_str[scales.index(scale)]
        subprocess.call(['mkdir','-p',dir_name])
        os.chdir(dir_name)
        for scale_diff in scales_diff:
            subprocess.call(['mkdir','-p',str(scale_diff)])
            os.chdir(str(scale_diff))
            
            for unit in units:
                subprocess.call(['mkdir','-p',str(unit)])
                os.chdir(str(unit))            
                subprocess.call(['scp','shurik@bugaboo.westgrid.ca:/home/shurik/IncreasingDiffusion/'+main_dir_name+"/"+str(scale)+"/"+str(scale_diff)+"/"+str(unit)+"/"+'density*.dat',"."])
                subprocess.call(['scp','shurik@bugaboo.westgrid.ca:/home/shurik/IncreasingDiffusion/'+main_dir_name+"/"+str(scale)+"/"+str(scale_diff)+"/"+str(unit)+"/"+"concentration.dat","."])
 
                os.chdir("..")
            os.chdir("..")
        os.chdir("..")
    os.chdir("..")

if __name__=="__main__":
    #modify_file()
    run_simulations("84")
