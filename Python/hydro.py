#!/usr/bin/python
import pylab
import numpy
import os
import math
def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]<=0 and prof[counter+1]>0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)
def Analyze_Air_Simulations():
    print os.getcwd()
    widths=[]
    velocities=[]
    re=[]
    
    for i in range(3,33,3):
        dir_name="AirResults/"+str(i)
        os.chdir(dir_name)
        
        name_phi="phase200000.dat"
        name_vel="velocity200000.dat"       
        phase=numpy.loadtxt(name_phi)
        vel=numpy.loadtxt(name_vel)        
        dims=phase.shape
        
        center=phase[dims[0]/2,:]
       
        z1 = numpy.min(numpy.where(center > 0.0))
        z2 = numpy.max(numpy.where(center > 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center<0.0))+dims[1]
            z1=numpy.max(numpy.where(center<0.0))
        print z1,z2
        print "Bubble length=",(z2-z1)/200.0
        print "Slug length=",15.0-(z2-z1)/200.0
        print "Liquid_velocity",numpy.sum(vel[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2)
        #print "Gas holdup",float(len(numpy.where(array['phi']<0)[0]))/(dims[1]*(dims[0]-2))
        prof=phase[:,((z1+z2)/2)%dims[1]]
        widths.append(Get_Zero(prof))     
        velocities.append(vel[dims[0]/2,z2%dims[1]])
        re.append(vel[dims[0]/2,z2%dims[1]]*dims[0]/(0.2/3.0))
    
      
        os.chdir("../..")
    
    capillaries=numpy.array(velocities)*(0.2/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    
    print "Widths=",widths
    print "Capillaries=",capillaries
    print "Velocities=",velocities
    
    
    #pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    #pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    fig1=pylab.figure()  
    pylab.plot(capillaries,widths,"go-",linewidth=3,markersize=10)
    fig2=pylab.figure()
    pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)
    pylab.xlim(0.02,1.5)
    pylab.ylim(ymin=0.01)
   
    #fig.subplots_adjust(left=0.15,bottom=0.15)  
    #pylab.xticks(fontsize=20)
    #pylab.yticks(fontsize=20)
    #pylab.xlabel(r'''$Ca$''',fontsize=30)
    #pylab.ylabel(r'''$\delta$''',fontsize=30)
    

if __name__=="__main__":
    Analyze_Simulations()
    pylab.show()