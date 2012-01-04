#!/usr/bin/python
import numpy
import pylab

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
    resources.vfYArray=numpy.linspace(0.0,1.0,len(ux[::2,1]))
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
 
     
    plot=Ngl.streamline(wks,uy[::2,::5],ux[::2,::5],resources)
    Ngl.contour(wks,geom,resources2)        
    #Ngl.contour(wks,phase,resources2)        
    #plot=Ngl.streamline(wks,vx,vy,resources)
    #Ngl.end()
    
    ### Producing outcome of velocities
    ux=numpy.loadtxt(dir_name+"vely0200000.dat")
    uy=numpy.loadtxt(dir_name+"velx0200000.dat")
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
    
    

if __name__=="__main__":
    file_dir="TweakingHydro2DRight/90/"
    print file_dir
    produce_streamlines(file_dir)
    #Ngl.end()
    
    
