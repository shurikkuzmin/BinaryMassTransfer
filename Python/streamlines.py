#!/usr/bin/python
import numpy
import pylab
import os
import math

def produce_bunch_streamlines(file_dir):
    for x in range(3,93,3):
        produce_streamlines(file_dir+str(x)+"/")
        print "Done with ",x

def produce_streamlines(dir_name):
    from PyNGL import Ngl
    geom=numpy.loadtxt(dir_name+"geometry.dat").transpose()
    ux=numpy.loadtxt(dir_name+"uy.dat").transpose()
    uy=numpy.loadtxt(dir_name+"ux.dat").transpose()
    dims=ux.shape
    #pylab.figure()
    #pylab.imshow(phase)
    #pylab.figure()
    #pylab.imshow(ux)
    #pylab.figure()
    #pylab.imshow(uy)    
    print dims
 
    y,x=numpy.mgrid[0:dims[0],0:dims[1]]        
        
    #positive=numpy.where(phase>0.0)
    #negative=numpy.where(phase<0.0)
 
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,dir_name+"original")
    resources = Ngl.Resources()
    
        
    #resources.tiMainFont    = "Times-Roman"
    #resources.tiMainOn=True
    #resources.tiMainString="Ca=0.22 "
           
    #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 3*0.25
    #resources.wkPaperSize="A5"
    #resources.nglFrame = False
    resources.vfXArray=numpy.linspace(0.0,15.0,len(ux[1,::5]))
    resources.vfYArray=numpy.linspace(0.0,1.0,len(ux[::5,1]))
    resources.nglFrame=False


    
    resources2=Ngl.Resources()
    resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources2.vpWidthF  = 3*0.25

    resources2.cnLineLabelsOn = False   # Turn off contour line labels.
        #resources2.cnLinesOn      = False   # Turn off contour lines.
    resources2.cnFillOn       = False    # Turn on contour fill.
    resources2.cnInfoLabelOn   = False 
  
        
    resources2.cnLevelSelectionMode = "ExplicitLevels"  # Select contour levels. 
    resources2.cnMinLevelValF       = 0.0
    #resources2.cnMaxLevelValF       = 0.001
    #resources2.cnLevelSpacingF      = 0.0
    resources2.cnLevelCount=1
    resources2.cnLevels=[0.0]
    #resources2.cnLineThicknesses=[3]
    resources2.cnMonoLineThickness=True
    resources2.cnLineThicknessF=3.0
    resources2.sfXArray=numpy.linspace(0.0,15.0,len(ux[1,:]))
    resources2.sfYArray=numpy.linspace(0.0,1.0,len(ux[:,1]))
 
     
    plot=Ngl.streamline(wks,uy[::5,::5],ux[::5,::5],resources)
    Ngl.contour(wks,geom,resources2)        
    #Ngl.contour(wks,phase,resources2)        
    #plot=Ngl.streamline(wks,vx,vy,resources)
    #Ngl.end()
    
    ### Producing outcome of velocities
    ux=numpy.loadtxt(dir_name+"velx0200000.dat").transpose()
    uy=numpy.loadtxt(dir_name+"vely0200000.dat").transpose()
    dims=ux.shape
    #pylab.figure()
    #pylab.imshow(phase)
    #pylab.figure()
    #pylab.imshow(ux)
    #pylab.figure()
    #pylab.imshow(uy)    
    print dims
 
    x,y=numpy.mgrid[0:dims[0],0:dims[1]]        
        
    #positive=numpy.where(phase>0.0)
    #negative=numpy.where(phase<0.0)
 
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,dir_name+"performed")
    ##resources = Ngl.Resources()
        
    #resources.tiMainFont    = "Times-Roman"
    #resources.tiMainOn=True
    #resources.tiMainString="Ca=0.22 "
           
    #resources.tiXAxisString = "streamlines"
    ##resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    ##resources.vpWidthF  = 3*0.25
    #resources.wkPaperSize="A5"
    ##resources.vfXArray=numpy.linspace(0.0,15.0,len(ux[1,::5]))
    ##resources.vfYArray=numpy.linspace(0.0,1.0,len(ux[::5,1]))
    ##resources.nglFrame = False

    
    plot=Ngl.streamline(wks,uy[::5,::5],ux[::5,::5],resources)
    Ngl.contour(wks,geom,resources2)        

    #Ngl.contour(wks,phase[::5,::50],resources2)        
    #Ngl.contour(wks,phase,resources2)        
    #plot=Ngl.streamline(wks,vx,vy,resources)

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)
    
def get_velocities_compare():

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


def get_velocities_compare_selection():

    print os.getcwd()
    widths=[]
    velocities=[]
    re=[]
    
    for i in [9,21,42,60,84]:
        dir_name="HydroResultsRight/"+str(i)
        os.chdir(dir_name)
        
        name_phi="phase300000.dat"
        name_velx="velocityx300000.dat"
        name_vely="velocityy300000.dat"
        
        phase=numpy.loadtxt(name_phi)
        velx=numpy.loadtxt(name_velx)
        vely=numpy.loadtxt(name_vely)

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





def produce_one_streamline():
    from PyNGL import Ngl

    dir_name="HydroResultsRight/"+str(3)
    os.chdir(dir_name)
    
    name_phi="phase300000.dat"
    name_velx="velocityx300000.dat"
    name_vely="velocityy300000.dat"
    
    phase=numpy.loadtxt(name_phi)
    velx=numpy.loadtxt(name_velx)
    vely=numpy.loadtxt(name_vely)

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
    print "Bubble velocity",velx[dims[0]/2,z2%dims[1]]        
    velx=velx-velx[dims[0]/2,z2%dims[1]]    
    #print "Gas holdup",float(len(numpy.where(array['phi']<0)[0]))/(dims[1]*(dims[0]-2))
    #prof=phase[:,((z1+z2)/2)%dims[1]]
    #widths.append(Get_Zero(prof))     
    #velocities.append()

    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,"initial")
    resources = Ngl.Resources()
    
        
    #resources.tiMainFont    = "Times-Roman"
    #resources.tiMainOn=True
    #resources.tiMainString="Ca=0.22 "
           
    #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 3*0.25
    #resources.wkPaperSize="A5"
    #resources.nglFrame = False
    resources.vfXArray=numpy.linspace(0.0,15.0,len(velx[1,::5]))
    resources.vfYArray=numpy.linspace(0.0,1.0,len(velx[::5,1]))
    resources.nglFrame=False


    
    resources2=Ngl.Resources()
    resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources2.vpWidthF  = 3*0.25

    resources2.cnLineLabelsOn = False   # Turn off contour line labels.
        #resources2.cnLinesOn      = False   # Turn off contour lines.
    resources2.cnFillOn       = False    # Turn on contour fill.
    resources2.cnInfoLabelOn   = False 
  
        
    resources2.cnLevelSelectionMode = "ExplicitLevels"  # Select contour levels. 
    resources2.cnMinLevelValF       = 0.0
    #resources2.cnMaxLevelValF       = 0.001
    #resources2.cnLevelSpacingF      = 0.0
    resources2.cnLevelCount=1
    resources2.cnLevels=[0.0]
    #resources2.cnLineThicknesses=[3]
    resources2.cnMonoLineThickness=True
    resources2.cnLineThicknessF=3.0
    resources2.sfXArray=numpy.linspace(0.0,15.0,len(velx[1,:]))
    resources2.sfYArray=numpy.linspace(0.0,1.0,len(velx[:,1]))
 
     
    plot=Ngl.streamline(wks,velx[::5,::5],vely[::5,::5],resources)
    Ngl.contour(wks,phase,resources2)        
    #Ngl.contour(wks,phase,resources2)        
    #plot=Ngl.streamline(wks,vx,vy,resources)
    #Ngl.end()
    
def comparison_bulk_velocities(file_dir):
    original_dir="TweakingHydro2DRight/"+str(file_dir)+"/"
    geometry=numpy.loadtxt(original_dir+"geometry.dat").transpose()
    ux_orig=numpy.loadtxt(original_dir+"ux.dat").transpose()
    uy_orig=numpy.loadtxt(original_dir+"uy.dat").transpose()
    
    ux_sim=numpy.loadtxt(original_dir+"vely0200000.dat").transpose()
    uy_sim=numpy.loadtxt(original_dir+"velx0200000.dat").transpose()
    
    dims_orig=ux_orig.shape
    dims_sim=ux_sim.shape
    print "Original=",dims_orig
    print "Simulations=",dims_sim
    
    pylab.figure()
    pylab.plot(ux_orig[dims_orig[0]/2,:])
    pylab.plot(ux_sim[dims_sim[0]/2,:])
    
    pylab.figure()
    pylab.plot(ux_orig[:,dims_orig[1]/2])
    pylab.plot(ux_sim[:,dims_sim[1]/2])
    
    pylab.figure()
    pylab.plot(ux_orig[:,0])
    pylab.plot(ux_sim[:,0])
    
    #Calculate some numbers as bulk velocities
    bulk_orig=numpy.sum(ux_orig[:,0])/dims_orig[0]
    bulk_sim=numpy.sum(ux_sim[:,0])/dims_sim[0]
    print "Bulk velocity original=",bulk_orig
    print "Bulk velocity simulations=",bulk_sim
    print "Simulations over original=",bulk_sim/bulk_orig
   
    profile=numpy.array(numpy.where(geometry[:,dims_orig[1]/2]>0))
    profile=profile.ravel()
    
    film_orig=numpy.sum(ux_orig[profile,dims_orig[1]/2])/len(profile)
    film_sim=numpy.sum(ux_sim[profile,dims_sim[1]/2])/len(profile)
    
    print "Bulk film velocity original=",film_orig
    print "Bulk film velocity simulations=",film_sim
    print "Simulations over original=",film_sim/film_orig

def produce_alexandra_image():
    from PyNGL import Ngl    
    file_dir="/home/shurik/Documents/People/Sasha/Data/"    
    phase=numpy.loadtxt(file_dir+"phi0015.dat")
    velx=numpy.loadtxt(file_dir+"velx0015.dat")
    vely=numpy.loadtxt(file_dir+"vely0015.dat")
    pylab.figure()    
    pylab.imshow(phase)
    pylab.figure()
    pylab.imshow(velx)
    pylab.figure()
    pylab.imshow(vely)
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
    print "Bubble velocity",velx[dims[0]/2,z2%dims[1]]        
    velx=velx-velx[dims[0]/2,z2%dims[1]]    
    #print "Gas holdup",float(len(numpy.where(array['phi']<0)[0]))/(dims[1]*(dims[0]-2))
    #prof=phase[:,((z1+z2)/2)%dims[1]]
    #widths.append(Get_Zero(prof))     
    #velocities.append()

    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,"alexandra")
    resources = Ngl.Resources()
    
        
    #resources.tiMainFont    = "Times-Roman"
    #resources.tiMainOn=True
    #resources.tiMainString="Ca=0.22 "
           
    #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 3*0.25
    #resources.wkPaperSize="A5"
    #resources.nglFrame = False
    resources.vfXArray=numpy.linspace(0.0,15.0,len(velx[1,::5]))
    resources.vfYArray=numpy.linspace(0.0,1.0,len(velx[::5,1]))
    resources.nglFrame=False


    
    resources2=Ngl.Resources()
    resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources2.vpWidthF  = 3*0.25

    resources2.cnLineLabelsOn = False   # Turn off contour line labels.
        #resources2.cnLinesOn      = False   # Turn off contour lines.
    resources2.cnFillOn       = False    # Turn on contour fill.
    resources2.cnInfoLabelOn   = False 
  
        
    resources2.cnLevelSelectionMode = "ExplicitLevels"  # Select contour levels. 
    resources2.cnMinLevelValF       = 0.0
    #resources2.cnMaxLevelValF       = 0.001
    #resources2.cnLevelSpacingF      = 0.0
    resources2.cnLevelCount=1
    resources2.cnLevels=[0.0]
    #resources2.cnLineThicknesses=[3]
    resources2.cnMonoLineThickness=True
    resources2.cnLineThicknessF=3.0
    resources2.sfXArray=numpy.linspace(0.0,15.0,len(velx[1,:]))
    resources2.sfYArray=numpy.linspace(0.0,1.0,len(velx[:,1]))
 
     
    plot=Ngl.streamline(wks,velx[::5,::5],vely[::5,::5],resources)
    Ngl.contour(wks,phase,resources2)        
    
    
def produce_one_image(file_dir):
    geom=numpy.loadtxt(file_dir+"geometry.dat")
    ux=numpy.loadtxt(file_dir+"vely0200000.dat")
    uy=numpy.loadtxt(file_dir+"velx0200000.dat")
    pylab.figure()    
    pylab.imshow(geom)
    pylab.figure()
    pylab.imshow(ux)
    pylab.figure()
    pylab.imshow(uy)
def produce_one_mass_image():
    mass=numpy.loadtxt("ScalingPeclet/density0020000.dat")
    pylab.figure()
    pylab.imshow(mass)
    pylab.colorbar()

def produce_new_image():
    from PyNGL import Ngl    
    file_dir="/home/shurik/Documents/Projects/binary_mass_transfer/Python/HydroAlexandra/9/"    
    phase=numpy.loadtxt(file_dir+"phase300000.dat")
    velx=numpy.loadtxt(file_dir+"velocityx300000.dat")
    vely=numpy.loadtxt(file_dir+"velocityy300000.dat")
    pylab.figure()    
    pylab.imshow(phase)
    pylab.figure()
    pylab.imshow(velx)
    pylab.figure()
    pylab.imshow(vely)
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
    print "Bubble velocity",velx[dims[0]/2,z2%dims[1]]        
    velx=velx-velx[dims[0]/2,z2%dims[1]]    
    #print "Gas holdup",float(len(numpy.where(array['phi']<0)[0]))/(dims[1]*(dims[0]-2))
    #prof=phase[:,((z1+z2)/2)%dims[1]]
    #widths.append(Get_Zero(prof))     
    #velocities.append()

    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,"new_force")
    resources = Ngl.Resources()
    
        
    #resources.tiMainFont    = "Times-Roman"
    #resources.tiMainOn=True
    #resources.tiMainString="Ca=0.22 "
           
    #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 3*0.25
    #resources.wkPaperSize="A5"
    #resources.nglFrame = False
    resources.vfXArray=numpy.linspace(0.0,15.0,len(velx[1,::5]))
    resources.vfYArray=numpy.linspace(0.0,1.0,len(velx[::5,1]))
    resources.nglFrame=False


    
    resources2=Ngl.Resources()
    resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources2.vpWidthF  = 3*0.25

    resources2.cnLineLabelsOn = False   # Turn off contour line labels.
        #resources2.cnLinesOn      = False   # Turn off contour lines.
    resources2.cnFillOn       = False    # Turn on contour fill.
    resources2.cnInfoLabelOn   = False 
  
        
    resources2.cnLevelSelectionMode = "ExplicitLevels"  # Select contour levels. 
    resources2.cnMinLevelValF       = 0.0
    #resources2.cnMaxLevelValF       = 0.001
    #resources2.cnLevelSpacingF      = 0.0
    resources2.cnLevelCount=1
    resources2.cnLevels=[0.0]
    #resources2.cnLineThicknesses=[3]
    resources2.cnMonoLineThickness=True
    resources2.cnLineThicknessF=3.0
    resources2.sfXArray=numpy.linspace(0.0,15.0,len(velx[1,:]))
    resources2.sfYArray=numpy.linspace(0.0,1.0,len(velx[:,1]))
 
     
    plot=Ngl.streamline(wks,velx[::5,::5],vely[::5,::5],resources)
    Ngl.contour(wks,phase,resources2)        



if __name__=="__main__":
    file_dir="TweakingHydro2DRight/9/"
    #file_dir="ScalingPeclet/9/"    
    #print file_dir
    produce_streamlines(file_dir)
    #get_velocities_compare()
    #get_velocities_compare_selection()
    #produce_one_streamline()
    #comparison_bulk_velocities(3)
    #produce_one_image(file_dir)    
    #produce_one_mass_image()    
    #produce_alexandra_image()    
    #produce_new_image()    
    pylab.show()
    #Ngl.end()
    
    
