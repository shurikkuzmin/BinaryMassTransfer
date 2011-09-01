#!/usr/bin/python
import numpy
import pylab


def CalculateDivergence(name_phase,name_velx,name_vely):
    phase=numpy.loadtxt(name_phase)
    velx=numpy.loadtxt(name_velx)
    vely=numpy.loadtxt(name_vely)
    
    grad_velx=numpy.gradient(velx)
    grad_vely=numpy.gradient(vely)
    
    divergence=grad_velx[0]+grad_vely[1]
    
    pylab.figure()
    pylab.imshow(divergence)
    
    pylab.figure()
    pylab.imshow(phase)
    pylab.figure()
    pylab.plot(divergence[100,:])
    pylab.figure()
    pylab.plot(divergence[:,200])

def CalculateDivergenceSailfish(name):
    arr=numpy.load(name)    
    phase=arr['phi']
    velx=arr['v'][0]
    vely=arr['v'][1]
    
    grad_velx=numpy.gradient(velx)
    grad_vely=numpy.gradient(vely)
    
    divergence=grad_velx[0]+grad_vely[1]
    
    pylab.figure()
    pylab.imshow(divergence)
    
    pylab.figure()
    pylab.imshow(phase)
    pylab.figure()
    pylab.plot(divergence[100,:])
    pylab.figure()
    pylab.plot(divergence[:,200])
    
    
if __name__=="__main__":
    name_phase="tmp/phase150000.dat"
    name_velx="tmp/velocityx150000.dat"
    name_vely="tmp/velocityy150000.dat"
    #CalculateDivergence(name_phase,name_velx,name_vely)
    name="SailfishData/20/capillary200000.npz"    
    CalculateDivergenceSailfish(name)    
    pylab.show()
    