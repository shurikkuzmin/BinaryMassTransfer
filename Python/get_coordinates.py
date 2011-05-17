#!/usr/bin/python
import pylab
import numpy
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
    print (zero-0.5)/(len(prof)-2)

def Get_Bubble(file_name):
    
    #fig=pylab.figure()
    array=numpy.load(file_name)
    dims=array['phi'].shape
    print "Dimensions=",dims
    center=array['phi'][dims[0]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[0]
        z1=numpy.max(numpy.where(center>0.0))
    print z1,z2

        
    prof=array['phi'][:,((z1+z2)/2)%dims[0]]
    prof_domain=prof[1:len(prof)-1]
    pylab.figure()
    pylab.plot(prof)
    #delta_x=1.0/float(dims[1]-2)
    #x_domain=numpy.arange(0.5*delta_x,1.0,delta_x)
    #x=numpy.arange(0.0,float(di))/ratio
        #pylab.imshow(array['phi'])
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
    #pylab.plot(x_domain,prof_domain,"r+",markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
    Get_Zero(prof)
    pylab.figure()
    pylab.imshow(array['phi'])
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

    #fig.subplots_adjust(left=0.15,bottom=0.15)  
    #pylab.xticks(fontsize=20)
    #pylab.yticks(fontsize=20)
    #pylab.xlabel(r'''$\delta$''',fontsize=30)
    #pylab.ylabel(r'''$\phi$''',fontsize=30)
    #labels=[r'''$H_{\mathrm{eff}}='''+str(value-2)+r'''$''' for value in ny]
    #pylab.legend(labels)
    #pylab.xlim(xmax=0.15)
    
    
if __name__=="__main__":
    file_name="Example/grid350000.npz"
    Get_Bubble(file_name)
    pylab.show()