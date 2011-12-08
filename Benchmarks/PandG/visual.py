import numpy
import pylab
import time

def visualize(file_name):
    phase=numpy.loadtxt(file_name)
    pylab.imshow(phase)
    
    pylab.colorbar()
    
if __name__=="__main__":
    for counter in range(0,10000,100):
        file_name="density00"+"0"*(5-len(str(counter)))+str(counter)+".dat"
        phase=numpy.loadtxt(file_name)
        fig=pylab.figure(1)
        c_levels=[0.5]
        #pylab.imshow(phase,norm=None,vmin=0.0,vmax=1.0,cle)
        cs=pylab.contour(phase,levels=c_levels)
        pylab.clabel(cs)
        #pylab.colorbar()
        fig.show()
        fig.clear()
        time.sleep(1)
        
    #visualize(file_name)
    pylab.show()