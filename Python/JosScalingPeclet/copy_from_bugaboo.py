import numpy
import os
import subprocess


def copy_files(dir_name,scales):
    os.chdir(dir_name)
    for scale in scales:
        subprocess.call(['mkdir','-p',scale+"Symmetry"])
        os.chdir(scale+"Symmetry")
        subprocess.call(['scp','shurik@bugaboo.westgrid.ca:/home/shurik/SymmetryScalingPeclet/'+dir_name+"/"+scale+"/concentration.dat","."])
        #subprocess.call(['scp','shurik@bugaboo.westgrid.ca:/home/shurik/SymmetryScalingPeclet/'+dir_name+"/"+scale+"/density*.dat","."])
        os.chdir("..")
    os.chdir("..")


if __name__=="__main__":
    copy_files("9",["1","2","4","6","10","15","20","40"])
    copy_files("21",["15","20"])

    #copy_files("84",10,5)
