#!/usr/bin/python
import numpy
import pylab
import os
import math

def get_zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)


def analyze_simulations(main_dir):
    print os.getcwd()
    widths=[]
    velocities=[]
    re=[]
    
    for i in range(3,93,3):
        dir_name=main_dir+"/"+str(i)
        os.chdir(dir_name)
        
        name_phi="phase300000.dat"
        name_velx="velocityx300000.dat"
        name_vely="velocityy300000.dat"
        phase=numpy.loadtxt(name_phi)
        velx=numpy.loadtxt(name_velx)
        vely=numpy.loadtxt(name_vely)
        #grad_velx=numpy.gradient(velx)
        #grad_vely=numpy.gradient(vely)
        #divergence=grad_velx[0]+grad_vely[1]
        
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
        widths.append(get_zero(prof))     
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
    numpy.savetxt(main_dir+"/capillaries.txt",zip(capillaries,widths))
    #fig.subplots_adjust(left=0.15,bottom=0.15)  
    #pylab.xticks(fontsize=20)
    #pylab.yticks(fontsize=20)
    #pylab.xlabel(r'''$Ca$''',fontsize=30)
    #pylab.ylabel(r'''$\delta$''',fontsize=30)


def get_bubble(dir_name):
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


def produce_bunch(main_dir):
    for counter in range(3,93,3):
        get_bubble(main_dir+"/"+str(counter)+"/")

def show_fields(main_dir):
    name_phi=main_dir+"/phase300000.dat"
    name_velx=main_dir+"/velocityx300000.dat"
    name_vely=main_dir+"/velocityy300000.dat"
    phase=numpy.loadtxt(name_phi)
    velx=numpy.loadtxt(name_velx)
    vely=numpy.loadtxt(name_vely)
    pylab.figure()    
    pylab.imshow(phase)
    pylab.colorbar()
    pylab.figure()    
    pylab.imshow(velx)
    pylab.colorbar()
    pylab.figure()    
    pylab.imshow(vely)
    pylab.colorbar()
        

def show_streamlines(dir_name):
    from PyNGL import Ngl

    phase=numpy.loadtxt(dir_name+"geometry.dat")
    ux=numpy.loadtxt(dir_name+"ux.dat")
    uy=numpy.loadtxt(dir_name+"uy.dat")
    dims=ux.shape
    pylab.figure()
    pylab.imshow(phase)
    pylab.figure()
    pylab.imshow(ux)
    pylab.figure()
    pylab.imshow(uy)    
    print dims
 
    x,y=numpy.mgrid[0:dims[0],0:dims[1]]        
        
    positive=numpy.where(phase>0.0)
    negative=numpy.where(phase<0.0)
    
    #large=numpy.where(numpy.logical_or(x<20,x>dims[0]-20))
        
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,"test")
    resources = Ngl.Resources()
        
    #resources.tiMainFont    = "Times-Roman"
    #resources.tiMainOn=True
    #resources.tiMainString="Ca=0.22 "
           
    #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 3*0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 0.25
    #resources.wkPaperSize="A5"
    #resources.nglFrame = False
    #resources.vfXArray=numpy.linspace(0.0,1.0,len(ux[1,::50]))
    #resources.vfYArray=numpy.linspace(0.0,15.0,len(ux[::5,1]))
    
    plot=Ngl.streamline(wks,uy[::5,::2],ux[::5,::2],resources)
    #Ngl.contour(wks,phase[::5,::50],resources2)        
    #Ngl.contour(wks,phase,resources2)        
    #plot=Ngl.streamline(wks,vx,vy,resources)
    Ngl.end()


if __name__=="__main__":
    #analyze_simulations("HydroResultsRightAcceleration")
    #analyze_simulations("HydroMichal")    
    #analyze_simulations("HydroGinzburg")
    #show_fields("HydroResultsRightAcceleration/12")
    #produce_bunch("HydroMichal")    
    #produce_bunch("HydroGinzburg")
    #show_streamlines("HydroResultsRightAcceleration/12/")    
    show_streamlines("HydroMichal/30/")
    #show_streamlines("HydroGinzburg/9/")
    pylab.show()