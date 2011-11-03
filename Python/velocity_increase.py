#!/usr/bin/python
import numpy
import pylab

def show_fields(file_name):
    density=numpy.loadtxt(file_name)
    pylab.figure(1)    
    pylab.imshow(density)
    pylab.figure(2)
    pylab.plot(density[100,:])
    
    print density[100,0],density[100,-1]

if __name__=="__main__":
    file_name="../VelocityCode/density0002000.dat"
    show_fields(file_name)
    pylab.show()