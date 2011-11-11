#!/usr/bin/python
import numpy
import pylab
import scipy.special

coeff=32
omega=1.0/(0.5+0.5/coeff)
radius=39
diffusion=1.0/3.0*(1.0/omega-0.5)
num_terms=20

def read_file(file_name):
    concentration=numpy.loadtxt(file_name)
    pylab.imshow(concentration)
    

def read_files(file_dir):
    styles=["kv","k^","ks","ko","kD"]    
    non_time=[]    
    legs=[]    
    for num,counter in enumerate(range(200*coeff,2000*coeff,400*coeff)):
        file_name=file_dir+"density"+(7-len(str(counter)))*"0"+str(counter)+".dat"
        concentration=numpy.loadtxt(file_name)
        dims=concentration.shape
        #pylab.figure()
        #pylab.imshow(concentration)
        x=numpy.arange(0,dims[1])        
        x_rad=numpy.linspace(0.0,1.0,radius+1)

        
        bessel=numpy.zeros_like(x)
        bessel_zeros=scipy.special.jn_zeros(0,num_terms)
        non_time.append(diffusion*counter/(radius*radius))
        for term in range(0,num_terms):
            bessel=bessel+2/(bessel_zeros[term]*scipy.special.jn(1,bessel_zeros[term]))\
                *numpy.exp(-bessel_zeros[term]*bessel_zeros[term]*diffusion*counter/(radius*radius))\
                *scipy.special.jn(0,bessel_zeros[term]*x/radius)
        bessel=1-bessel
        #bessel[numpy.where(bessel)>1]=0
        #bessel=numpy.roll(bessel,dims[1]/2)
        print "Shape",dims        
        pylab.figure(1)
        pylab.plot(x_rad,concentration[(dims[0]-1)/2,(dims[1]-1)/2:(dims[1]-1)/2+radius+1],styles[num])        
        #pylab.plot(concentration[dims[0]/2-1,dims[1]/2:dims[1]/2+radius])
        #pylab.figure(2)
        pylab.plot(x_rad,bessel[numpy.where(x<=radius)],styles[num]+"--",markerfacecolor="None")
    legs=[[r'''$\tau='''+str(time)[0:6]+'''$''',r'''$\tau='''+str(time)[0:6]+'''$'''] for time in non_time]    
    legs=numpy.ravel(legs)   
    pylab.title(r'''$D='''+str(diffusion)[0:6]+'''$''',fontsize=30)
    pylab.legend(legs,fancybox=True,loc=2)
    pylab.xlabel(r'''$x$''',fontsize=20)
    pylab.ylabel(r'''$C$''',fontsize=20)
    #pylab.legend(legs,fancybox=True,labelspacing=0.1)
    pylab.savefig("cylinder"+str(diffusion)[2:6]+".eps",dpi=300)        
        
def show_bessel():
    x=numpy.arange(0.0,10.0,0.1)
    bessel=scipy.special.jn(0,x)
    pylab.plot(x,bessel)

def read_film(file_dir):
    concentration=numpy.loadtxt(file_dir+"conc_initial.dat")
    ux=numpy.loadtxt(file_dir+"ux_initial.dat")
    uy=numpy.loadtxt(file_dir+"uy_initial.dat")
    print concentration.shape
    pylab.figure()    
    pylab.imshow(concentration)
    pylab.figure()
    pylab.plot(concentration[:,200])
    pylab.figure()
    pylab.plot(concentration[:,0])
    pylab.figure()
    pylab.imshow(uy)
    pylab.colorbar()

if __name__=="__main__":
    #file_name="../Benchmarks/density0001000.dat"
    file_name="../Benchmarks/conc_initial.dat"    
    file_dir="../Benchmarks/"
    #read_file(file_name)
    
    #read_files(file_dir)
    #show_bessel()    
    read_film(file_dir)    
    pylab.show()
