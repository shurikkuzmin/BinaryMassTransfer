#!/usr/bin/python
import pylab
import numpy
import os
import math

#physical parameters
gamma=0.0728
diam=1.5e-3

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)
def Analyze_Simulations():
    print os.getcwd()
    widths=[]
    velocities=[]
    re=[]
    
    for i in range(3,93,3):
        dir_name="HydroResultsRight/"+str(i)
        os.chdir(dir_name)
        
        name_phi="phase300000.dat"
        name_velx="velocityx300000.dat"
        name_vely="velocityy300000.dat"
        phase=numpy.loadtxt(name_phi)
        velx=numpy.loadtxt(name_velx)
        vely=numpy.loadtxt(name_vely)
        grad_velx=numpy.gradient(velx)
        grad_vely=numpy.gradient(vely)
        divergence=grad_velx[0]+grad_vely[1]
        
        if numpy.max(divergence)>0.0001:
            pylab.figure()            
            pylab.imshow(divergence)
        dims=phase.shape
        
        center=phase[dims[0]/2,:]
       
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+dims[1]
            z1=numpy.max(numpy.where(center>0.0))
        print z1,z2
        print "Bubble length=",(z2-z1)/200.0
        print "Slug length=",15.0-(z2-z1)/200.0
        print "Liquid_velocity",numpy.sum(velx[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2)
        #print "Gas holdup",float(len(numpy.where(array['phi']<0)[0]))/(dims[1]*(dims[0]-2))
        prof=phase[:,((z1+z2)/2)%dims[1]]
        widths.append(Get_Zero(prof))     
        velocities.append(velx[dims[0]/2,z2%dims[1]])
        re.append(velx[dims[0]/2,z2%dims[1]]*dims[0]/(2.0/3.0))
    
      
        os.chdir("../..")
    
    capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    
    print "Widths=",widths
    print "Capillaries=",capillaries
    print "Velocities=",velocities
    
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

    
def Get_Bubble(dir_name):
    dir_name="HydroResults/"+dir_name+"/"
    phase=numpy.loadtxt(dir_name+"phase300000.dat")
    velx=numpy.loadtxt(dir_name+"velocityx300000.dat")
    vely=numpy.loadtxt(dir_name+"velocityy300000.dat")
    dims=phase.shape

    print "Dimensions=",dims

    center=phase[dims[0]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    slug_length=dims[1]-z2+z1    
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[1]
        z1=numpy.max(numpy.where(center>0.0))
        slug_length=z1-z2%dims[1]
    print z1,z2
    
    ##Performing perturbations
    
    ##reference frame
    interface_velocity=velx[dims[0]/2,z2%dims[1]]
    print "Interface velocity=",interface_velocity
    print "Capillary=",interface_velocity*2.0/3.0/math.sqrt(8*0.04*0.04/9.0)
    bubble_reference=velx-interface_velocity
 
    ##Rolling the bubble to the end with certain shift
    shift=slug_length/2
    print dims[1]-(z2%dims[1])
    bubble_reference=numpy.roll(bubble_reference,dims[1]-(z2%dims[1])-shift,axis=1)
    phase=numpy.roll(phase,dims[1]-(z2%dims[1])-shift,axis=1)    
    vely=numpy.roll(vely,dims[1]-(z2%dims[1])-shift,axis=1)    
 
    #pylab.figure()
    #pylab.plot(velx[:,((z2-z1)/2)%dims[1]])
    #pylab.figure()
    #pylab.plot(bubble_reference[:,dims[1]/2])
    ##pylab.figure()
    ##pylab.imshow(phase)
    
    #Fliping array
    bubble_reference=-numpy.fliplr(bubble_reference)
    phase=numpy.fliplr(phase)
    vely=numpy.fliplr(vely)
    #Calculate_Perimeter(phase)
    ##pylab.figure()
    ##pylab.imshow(bubble_reference)
    
    ##pylab.figure()
    ##pylab.plot(bubble_reference[(dims[0]-1)/2,:])
    
    ##pylab.figure()
    ##pylab.plot(bubble_reference[:,0])
    
    
    positive=numpy.where(phase>0.0)
    negative=numpy.where(phase<=0.0)
    geometry=numpy.zeros([dims[0],dims[1]],dtype="int")    
    geometry[positive]=1
    geometry[negative]=-1
    
    cx=[0,1,0,-1,0,1,-1,-1,1]
    cy=[0,0,1,0,-1,1,1,-1,-1]
    dirs=zip(cx,cy)
    for x,y in zip(negative[0],negative[1]):
        for dirx,diry in dirs:
            posx=x+dirx
            posy=y+diry
            if geometry[posx,posy]==1:
                geometry[x,y]=0
                break
            
    #pylab.figure()
    #pylab.imshow(geometry)
    
    ##numpy.savetxt("geometry.dat",geometry,fmt="%d")
    ##numpy.savetxt("ux.dat",bubble_reference)
    ##numpy.savetxt("uy.dat",vely)
    numpy.savetxt(dir_name+"/geometry.dat",geometry.transpose(),fmt="%d")
    numpy.savetxt(dir_name+"/ux.dat",bubble_reference.transpose())
    numpy.savetxt(dir_name+"/uy.dat",vely.transpose())

def Produce_Bunch():
    for counter in range(3,93,3):
        Get_Bubble(str(counter))
        
      
        
def Produce_Mass_Bunch():
    print os.getcwd()
 
    styles=['kv','ks','ko','k^','k>','k<','kD','kp']
    aver_coefficient=[]
    bubble_velocities=[]
    
    fig=pylab.figure(99)

    for counter in range(3,39,3):    
        dir_name="HydroResults/"+str(counter)+"/"
        phase=numpy.loadtxt(dir_name+"phase300000.dat")
        velx=numpy.loadtxt(dir_name+"velocityx300000.dat")
        vely=numpy.loadtxt(dir_name+"velocityy300000.dat")
        dims=phase.shape

        print "Dimensions=",dims

        center=phase[dims[0]/2,:]
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        slug_length=dims[1]-z2+z1    
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+dims[1]
            z1=numpy.max(numpy.where(center>0.0))
            slug_length=z1-z2%dims[1]
        print z1,z2

        interface_velocity=velx[dims[0]/2,z2%dims[1]]
        capillary=interface_velocity*(2.0/3.0)/math.sqrt(8*0.04*0.04/9.0)
        reynolds=interface_velocity*dims[0]/(2.0/3.0)        
        phys_velocity=math.sqrt(capillary*reynolds*gamma/(1000*diam))        
        
        bubble_velocities.append(phys_velocity)
        
        print "Interface velocity=",interface_velocity
        print "Capillary=",capillary
        print "Reynolds=",reynolds
        
        dir_name="HydroResults/Mass/"+str(counter)+"/"        
        name=dir_name+"concentration.dat"
        array=numpy.loadtxt(name)

        #print array
        new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        #print "Data=",new_conc
        deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)
        
        coefficient=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        pylab.figure(99)
        pylab.plot(array[1:,0],coefficient)#,styles[counter])
        aver_coefficient.append(numpy.mean(coefficient[len(coefficient)/2:]))
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)
    #pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    print "Aver_coefficient=",aver_coefficient
    return aver_coefficient



if __name__=="__main__":
    Analyze_Simulations()
    #Produce_Hydro_Bunch()    
    #Produce_Mass_Bunch()
    pylab.show()