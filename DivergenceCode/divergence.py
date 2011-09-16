#!/usr/bin/python
import numpy
import pylab

if __name__=="__main__":
    arr=numpy.array([[1,10],[3,18],[5,12]])
    print arr.shape    
    grad=numpy.gradient(arr)
    print grad
    
    velx=numpy.loadtxt("ux.dat")
    vely=numpy.loadtxt("uy.dat")
    gradx=numpy.gradient(velx)
    grady=numpy.gradient(vely)
    pylab.figure()    
    pylab.imshow(velx)
    pylab.figure()
    pylab.imshow(vely)
    
    pylab.figure()
    pylab.imshow(gradx[0]+grady[1])
    #pylab.figure()
    #pylab.imshow(gradx[1]+grady[0])
    print velx.shape
    pylab.show()
    pylab.savetxt("divergence.dat",gradx[0]+grady[1])