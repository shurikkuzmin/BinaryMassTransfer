#import numpy
import os
import subprocess

def copy_files():
    capillary_str=[str(x) for x in range(9,93,10)]
    for length in range(200,1400,150):
        subprocess.call(['mkdir','-p',str(length)])
        os.chdir(str(length))
        for i in range(0, len(capillary_str)):
            dir_temp=capillary_str[i]
            subprocess.call(['mkdir','-p',dir_temp])
            os.chdir(dir_temp)
            subprocess.call(["scp","shurik@bugaboo.westgrid.ca:/home/shurik/ShowCases/"+str(length)+"/"+capillary_str[i]+"/"+"tmp/*300000.dat","."])
            os.chdir("..")
        os.chdir("..")

if __name__=="__main__":
    copy_files()
