import numpy
import os
import subprocess

def copy_files(server):
    subprocess.call(['scp','binary.cpp',server+".westgrid.ca:/home/shurik/ShowCasesScaled/"]) 
    subprocess.call(['scp','binary.pbs',server+".westgrid.ca:/home/shurik/ShowCasesScaled/"]) 
    subprocess.call(['scp','govern.py',server+".westgrid.ca:/home/shurik/ShowCasesScaled/"]) 

if __name__=="__main__":
    copy_files("bugaboo")
