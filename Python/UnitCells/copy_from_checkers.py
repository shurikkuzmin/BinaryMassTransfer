import numpy
import os
import subprocess

def run_simulations():
    capillary_str=[str(x) for x in [9,21]]
    scale=[2,4,8,10,20,40]
    scale_str=[str(x) for x in [2,4,8,10,20,40]]
    for counter,dir_temp in enumerate(capillary_str):
        os.chdir(dir_temp)
        for scale_counter in range(0,len(scale_str)-counter):
            #subprocess.call(['mkdir','-p',scale_str[scale]])
           
            os.chdir(scale_str[scale_counter])
            name=str(20000*scale[len(scale_str)-counter-1]/scale[scale_counter])
            print name
            #subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/ScalingPeclet/'+dir_temp+"/"+scale_str[scale]+"/concentration.dat","."])
            subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/ScalingPeclet/'+dir_temp+"/"+scale_str[scale_counter]+"/density"+'0'*(7-len(name))+name+".dat","."])
            os.chdir("..")
        os.chdir("..")


def copy_files(dir_name,factor,subtract):
    scale=[2,4,8,10,20,40]
    scale_str=[str(x) for x in [2,4,8,10,20,40]]
    os.chdir(dir_name)
    for scale_counter in range(0,len(scale_str)-subtract):
        #subprocess.call(['mkdir','-p',scale_str[scale]])

        os.chdir(scale_str[scale_counter])
        name=str(factor*20000*scale[len(scale_str)-subtract-1]/scale[scale_counter])
        print name
            #subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/ScalingPeclet/'+dir_temp+"/"+scale_str[scale]+"/concentration.dat","."])
        subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/ScalingPeclet/'+dir_name+"/"+scale_str[scale_counter]+"/density"+'0'*(7-len(name))+name+".dat","."])
        os.chdir("..")
    os.chdir("..")


if __name__=="__main__":
    #modify_file()
    #run_simulations()
    #copy_files("42",5,4)
    #copy_files("60",10,5)
    copy_files("84",10,5)
