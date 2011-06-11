#!/usr/bin/python
import pylab
import numpy
import math
import os

def Streamlines_Plot(ux,uy,phase):
    from PyNGL import Ngl
    from numpy import genfromtxt

    dims=phase.shape        
    x,y=numpy.mgrid[0:dims[0],0:dims[1]]        
  
    center=phase[(dims[0]-1)/2,:]
           
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[0]
        z1=numpy.max(numpy.where(center>0.0))
    print z1,z2
    
    print "Capillary=",ux[(dims[0]-1)/2,z2%dims[1]]*2.0/3.0/math.sqrt(8.0*0.04*0.04/9.0)
    positive=numpy.where(phase>0.0)
    negative=numpy.where(phase<0.0)
    
    large=numpy.where(numpy.logical_or(x<20,x>dims[0]-20))
    ux=ux-ux[(dims[0]-1)/2,z2%dims[1]]        
    #ux[negative]=None
    #ux[large]=None
    #vz_diff_mask[bounds]=None
        
    ux_mask=ux[::5,::100]
    uy_mask=uy[::5,::100]
    x_mask=x[::5,::100]
    y_mask=y[::5,::100]

        
    #pylab.figure(figsize=(20,10))
    #pylab.quiver(y_mask,x_mask,ux_mask,uy_mask,headwidth=6,minlength=0.1)
    #pylab.contour(array['phi'],[0.0],linewidths=[3])
        
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,"vortex")
    resources = Ngl.Resources()
     
     
        
    #resources.tiMainFont    = "Times-Roman"
    #resources.tiMainOn=True
    #resources.tiMainString="Ca=0.22 "
           
    #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 3*0.25
    resources.wkPaperSize="A5"
    resources.nglFrame = False
    resources.vfXArray=numpy.linspace(0.0,15.0,len(ux[1,::50]))
    resources.vfYArray=numpy.linspace(0.0,1.0,len(ux[::5,1]))


        
        
    resources2=Ngl.Resources()
    #resources2.tiMainFont    = "Times-Roman"
    #resources2.tiXAxisString = "streamlines"
    #resources2.tiMainOn=True
    #resources2.tiMainString="Ca=0.22"
        
    resources2.wkPaperSize="A5"
    resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources2.vpWidthF  = 3*0.25
    resources2.nglFrame = False
        
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
        
    resources2.lbLabelBarOn=False
    resources2.lbLabelsOn=False
    resources2.sfXArray=numpy.linspace(0.0,15.0,len(ux[1,:]))
    resources2.sfYArray=numpy.linspace(0.0,1.0,len(ux[:,1]))
        
    plot=Ngl.streamline(wks,ux[::5,::50],uy[::5,::50],resources)
    #Ngl.contour(wks,phase[::5,::50],resources2)        
    Ngl.contour(wks,phase,resources2)        
    #plot=Ngl.streamline(wks,vx,vy,resources)
    Ngl.end()

def Contour_Extraction():
    print "Nothing"   
   #for counter,omega_value in enumerate(omega_int):
   #     pylab.figure(1) 
   #     print "Omega=",omega_value 
   #     filename="PhaseTRTScript/domainomega"+omega_str[omega.index(omega_value)]+"gamma"+gamma_str+"ux"+ux_str[ux.index(ux_value)]+".dat"
   #     print "Filename=",filename        
   #     data=numpy.loadtxt(filename)
   #     #pylab.figure()        
   #     CS=pylab.contour(x,y,data,[0.5]) #colors=[colors[counter]],linewidths=[4],linestyles=styles[counter])
   #     path0=CS.collections[0].get_paths()[0]
   #     cont=path0.vertices
   #     pylab.figure(2)
   #     #if counter<3:
   #     pylab.plot(cont[:,0],cont[:,1],styles[counter],markerfacecolor="None",markersize=counter+10,alpha=0.5)

   #probe=vtk.vtkProbeFilter()
   # probe.SetInputConnection(contour.GetOutputPort())    
   # probe.SetSource(image)
   # probe.SpatialMatchOn()    
   # probe.Update()

   # print "Probe=",probe.GetOutput()

   # cont=probe.GetOutput()
   # vel=cont.GetPointData().GetArray("Velocity")    
   # phi=cont.GetPointData().GetArray("Phase")    
   # cont_points=cont.GetPoints()
   # x_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
   # y_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
   # z_numpy=numpy.zeros(cont_points.GetNumberOfPoints())    
    
   # velx_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
   # vely_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
   # velz_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
    
   # phi_numpy=numpy.zeros(cont_points.GetNumberOfPoints())

   # for counter in range(0,cont.GetPoints().GetNumberOfPoints()):
   #     x,y,z=cont_points.GetPoint(counter)
   #     x_numpy[counter]=x
   #     y_numpy[counter]=y
   #     z_numpy[counter]=z
   #     velx_numpy[counter]=vel.GetTuple3(counter)[0]
   #     vely_numpy[counter]=vel.GetTuple3(counter)[1]
   #     velz_numpy[counter]=vel.GetTuple3(counter)[2]
   #     phi_numpy[counter]=phi.GetTuple1(counter) 

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)

def Calculate_Perimeter(phase):
       #for counter,omega_value in enumerate(omega_int):
   #     pylab.figure(1) 
   #     print "Omega=",omega_value 
   #     filename="PhaseTRTScript/domainomega"+omega_str[omega.index(omega_value)]+"gamma"+gamma_str+"ux"+ux_str[ux.index(ux_value)]+".dat"
   #     print "Filename=",filename        
   #     data=numpy.loadtxt(filename)
   #     #pylab.figure()        
   #     CS=pylab.contour(x,y,data,[0.5]) #colors=[colors[counter]],linewidths=[4],linestyles=styles[counter])
   #     path0=CS.collections[0].get_paths()[0]
   #     cont=path0.vertices
   #     pylab.figure(2)
   #     #if counter<3:
   #     pylab.plot(cont[:,0],cont[:,1],styles[counter],markerfacecolor="None",markersize=counter+10,alpha=0.5)
    
    dims=phase.shape
    #x,y=numpy.mgrid[0:dims[0],0:dims[1]]
    pylab.figure(99)
    CS=pylab.contour(phase,[0.0])
    path0=CS.collections[0].get_paths()[0]
    cont=path0.vertices
    perimeter=0    
    for i,(x,y) in enumerate(cont[:-1]):
        perimeter+=math.sqrt((x-cont[i+1,0])**2+(y-cont[i+1,1])**2)
    print "Perimeter=",perimeter
    print "a=",perimeter/((dims[0]-2)*dims[1]*1.5e-3);
    #pylab.plot(cont[:,0],cont[:,1],"+")


def Get_Bubble(file_name,dir_name):
    
    #fig=pylab.figure()
    array=numpy.load(file_name)
    dims=array['phi'].shape
    velx=array['v'][0]
    vely=array['v'][1]
    phase=array['phi']
    print "Dimensions=",dims
    print array.files
    center=phase[(dims[0]-1)/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[1]
        z1=numpy.max(numpy.where(center>0.0))
    print z1,z2

        
    
    #Performing perturbations
    #reference frame
    interface_velocity=velx[(dims[0]-1)/2,z2%dims[1]]
    print "Interface velocity=",interface_velocity
    print "Capillary=",interface_velocity*2.0/3.0/math.sqrt(8*0.04*0.04/9.0)
    bubble_reference=array['v'][0]-interface_velocity
 
    #Streamlines_Plot(velx,vely,phase)

    
    #Rolling the bubble to the end with certain shift
    shift=50
    bubble_reference=numpy.roll(bubble_reference,dims[1]-(z2%dims[1])-shift)
    phase=numpy.roll(phase,dims[1]-(z2%dims[1])-shift)    
    vely=numpy.roll(vely,dims[1]-(z2%dims[1])-shift)    
    #pylab.figure()
    #pylab.plot(bubble_reference[:,z2%dims[1]+30])
    #pylab.figure()
    #pylab.plot(bubble_reference[(dims[0]-1)/2,:])
    #pylab.figure()
    #pylab.imshow(phase)
    
    #Fliping array
    bubble_reference=-numpy.fliplr(bubble_reference)
    phase=numpy.fliplr(phase)
    vely=numpy.fliplr(vely)
    Calculate_Perimeter(phase)
    #pylab.figure()
    #pylab.imshow(bubble_reference)
    
    #pylab.figure()
    #pylab.plot(bubble_reference[(dims[0]-1)/2,:])
    
    #pylab.figure()
    #pylab.plot(bubble_reference[:,0])
    
    
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
            
    pylab.figure()
    pylab.imshow(geometry)
    
    #numpy.savetxt("geometry.dat",geometry,fmt="%d")
    #numpy.savetxt("ux.dat",bubble_reference)
    #numpy.savetxt("uy.dat",vely)
    numpy.savetxt(dir_name+"/geometry.dat",geometry.transpose(),fmt="%d")
    numpy.savetxt(dir_name+"/ux.dat",bubble_reference.transpose())
    numpy.savetxt(dir_name+"/uy.dat",vely.transpose())

def Produce_Bunch():
    print os.getcwd()
    capillary_str=["3","5","8","10","20","40","60","80"]

    velocities_real=[]
    re_real=[]
    widths_real=[]

    for dir_temp in capillary_str:
        dir_name="../../binary_microchannel/Sailfish/Capillary/Results/"+dir_temp
        name=dir_name+"/capillary200000.npz"
        if not os.path.isdir(dir_temp):
            os.mkdir(dir_temp)
        Get_Bubble(name,dir_temp)
        array=numpy.load(name)
        #prof=array['phi'][:,exam[i]]
        dims=array['phi'].shape
        print dims
        #x=numpy.arange(0.0,float(ny[i]))/ratio
        #pylab.imshow(array['phi'])
    
        #pylab.plot(prof)
        #widths.append(Get_Zero(prof))
        vel=array['v'][0]
        center=array['phi'][dims[0]/2,:]
       
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+dims[1]
            z1=numpy.max(numpy.where(center>0.0))
        print z1,z2
        
        prof_real=array['phi'][:,((z1+z2)/2)%dims[1]]
        widths_real.append(Get_Zero(prof_real))     
        
        
        #vel_prof=vel[:,exam[i]]
        #velocities.append(vel_prof[len(vel_prof)/2])
        velocities_real.append(vel[dims[0]/2,z2%dims[1]])
        re_real.append(vel[dims[0]/2,z2%dims[1]]*dims[0]/(2.0/3.0))
    
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
        #pylab.plot(x,prof,style_diff[i],markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        #Get_Zero(prof)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        #os.chdir("../..")
    
    fig=pylab.figure()
    #capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    capillaries_real=numpy.array(velocities_real)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    
    
    
if __name__=="__main__":
    file_name="Example/capillary200000.npz"
    #Get_Bubble(file_name)
    Produce_Bunch()    
    pylab.show()
