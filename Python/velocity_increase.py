#!/usr/bin/python
import numpy
import pylab

def show_fields(dir_name,file_name):
    density=numpy.loadtxt(dir_name+file_name)
    ux=numpy.loadtxt(dir_name+"vely0080000.dat").transpose()    
    pylab.figure(1)    
    pylab.imshow(density)
    pylab.colorbar()
    pylab.figure(2)
    pylab.plot(density[100,:])
    
    abstract_array=numpy.zeros_like(density)
    abstract_array[numpy.where(density>1.1)]=1
    pylab.figure(3)
    pylab.imshow(abstract_array)
    pylab.figure(4)
    pylab.imshow(ux)
    pylab.colorbar()
    
    
    print density[100,0],density[100,-1]
    print density[numpy.where(density>1.1)]
    print numpy.where(density>1.1)

#To be transfered to C++ file
def prepare_files(dir_name):
    geometry=numpy.loadtxt(dir_name+"geometry.dat").transpose()    
    ux=numpy.loadtxt(dir_name+"ux.dat").transpose()
    uy=numpy.loadtxt(dir_name+"uy.dat").transpose()
    dims=geometry.shape

    pylab.figure(1)
    pylab.imshow(geometry)
    pylab.figure(2)
    pylab.imshow(ux)
    print numpy.where(ux<-0.01)
    pylab.colorbar()
    pylab.title("Original")
    pylab.figure(3)
    pylab.plot(ux[:,100])
    pylab.figure(4)
    pylab.plot(ux[dims[0]/2,:])
    
    pylab.figure(99)
    pylab.imshow(uy)
    pylab.colorbar()
    pylab.title("Original")
    pylab.figure(199)
    pylab.plot(uy[dims[0]/2,:])

    interface_velocity=ux[0,0]

   
    negative=numpy.where(geometry==-1)
    zerolevel=numpy.where(geometry==0)
    ux[negative]=0.0
    uy[negative]=0.0
    ux[zerolevel]=0.0
    #ux[0,:]=interface_velocity
    #ux[-1,:]=interface_velocity
    uy[zerolevel]=0.0
    
    pylab.figure(6)
    pylab.imshow(ux)
    pylab.colorbar()
    pylab.title("PostProcessed")
    pylab.figure(100)
    pylab.imshow(uy)
    pylab.colorbar()
    pylab.title("PostProcessed")
    pylab.figure(199)
    pylab.plot(uy[dims[0]/2,:])
    
    pylab.figure(4)
    pylab.plot(ux[dims[0]/2,:])
    ux=ux.transpose()
    uy=uy.transpose()
    geometry=geometry.transpose()

    
    #returning back geometry
    pylab.savetxt(dir_name+"ux_non.dat",ux)
    pylab.savetxt(dir_name+"uy_non.dat",uy)
    pylab.savetxt(dir_name+"geometry_non.dat",geometry,fmt="%d")
    
    
def compare_hydro_fields(file_dir,file_name):
    ux_orig=numpy.loadtxt(file_dir+"ux_non.dat").transpose()
    ux_proc=numpy.loadtxt(file_dir+file_name).transpose()    
    dims=ux_orig.shape    
    print "Compare", dims, ux_proc.shape    
    pylab.figure(1)    
    pylab.imshow(ux_orig)
    pylab.colorbar()
    pylab.title("Original")    
    pylab.figure(2)
    pylab.imshow(ux_proc)
    pylab.colorbar()
    pylab.title("Processed")
    pylab.figure(3)
    pylab.plot(ux_orig[dims[0]/2,:])
    #pylab.title("Original")
    #pylab.figure(4)
    pylab.plot(ux_proc[dims[0]/2,:])
    #pylab.title("Processed")
    pylab.figure(5)
    pylab.plot(ux_orig[:,100])
    pylab.plot(ux_proc[:,100])
    
    
if __name__=="__main__":
    dir_name="../VelocityCode/"
    file_name="density0003200.dat"
    #file_name="SailfishMassCentralPeriodic/80/density3960000.dat"    
    show_fields(dir_name,file_name)
    #file_dir="../VelocityCode/"
    #file_name="vely0080000.dat"    
    #compare_hydro_fields(file_dir,file_name)    
    #prepare_files(file_dir)    
    pylab.show()