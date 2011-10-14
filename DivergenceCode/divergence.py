#!/usr/bin/python
import numpy
import pylab
from numpy import ma

def poisson(div):
    dims=div.shape
    phi=numpy.zeros(dims)
    
    for counter in range(1,100):
        phi[:,0]=phi[:,1]
        phi[:,-1]=phi[:,-2]
        for i in range(0,dims[0]):
            for j in range(1,dims[1]-1):
                topi=(i+1+dims[0])%dims[0]
                bottomi=(i-1+dims[0])%dims[0]
                delta=0.25*(phi[topi,j]+phi[bottomi,j]+phi[i,j+1]+phi[i,j-1]-4*phi[i,j]+4*div[i,j])
                phi[i,j]=phi[i,j]+delta
        print "Counter=",counter
    return numpy.gradient(phi)

def read_poisson(velx,vely):
    divx=numpy.loadtxt("divx.dat")
    divy=numpy.loadtxt("divy.dat")
    resultx=numpy.gradient(velx-divx)
    resulty=numpy.gradient(vely-divy)
    pylab.figure()
    pylab.imshow(resultx[0]+resulty[1])
    pylab.colorbar()
def fourier(velx,vely):
    import numpy.fft
    #make it symmetric
    dims=velx.shape
    
    velx[dims[0]/2:dims[0]-1,:]=numpy.fliplr(velx)[dims[0]/2:dims[0]-1,:]
    vely[dims[0]/2:dims[0]-1,:]=numpy.fliplr(vely)[dims[0]/2:dims[0]-1,:]
    fftx=numpy.fft.fft2(velx)
    ffty=numpy.fft.fft2(vely)
    
    resultx=numpy.zeros(velx.shape,dtype="complex")
    resulty=numpy.zeros(vely.shape,dtype="complex")
    resultx[0,0]=fftx[0,0]
    resulty[0,0]=ffty[0,0]
    for coorx in range(0,dims[0]):
        for coory in range(0,dims[1]):
            if coorx==0 and coory==0:
                continue
            freqx=float(coorx)/float(dims[0])
            freqy=float(coory)/float(dims[1])
            resultx[coorx,coory]=fftx[coorx,coory]-(fftx[coorx,coory]*freqx+ffty[coorx,coory]*freqy)*freqx/(freqx*freqx+freqy*freqy)
            resulty[coorx,coory]=ffty[coorx,coory]-(fftx[coorx,coory]*freqx+ffty[coorx,coory]*freqy)*freqy/(freqx*freqx+freqy*freqy)
            tangentialx[coorx,coory]=(fftx[coorx,coory]*freqx+ffty[coorx,coory]*freqy)*freqx/(freqx*freqx+freqy*freqy)
            tangentialy[coorx,coory]=(fftx[coorx,coory]*freqx+ffty[coorx,coory]*freqy)*freqy/(freqx*freqx+freqy*freqy)
            
    #print fft[0,0],fft[2999,0]
    #antifft=numpy.fft.ifft2(fftx)
    antidivx=numpy.fft.ifft2(resultx)
    antidivy=numpy.fft.ifft2(resulty)
    print "Fouriered velx=",antidivx
    print "Fouriered vely=",antidivy
    #print antifft.shape
    #print numpy.real(antifft)
    #pylab.figure()
    #pylab.imshow(numpy.real(antifft))
    #pylab.colorbar()
    pylab.figure()
    pylab.imshow(numpy.real(antidivx))
    pylab.colorbar()
    pylab.figure()
    pylab.imshow(numpy.real(antidivy))
    pylab.colorbar()
    pylab.savetxt("ux_non.dat",numpy.real(antidivx))
    pylab.savetxt("uy_non.dat",numpy.real(antidivy))
    return numpy.real(antidivx),numpy.real(antidivy)    

def calculate_center_fields():
    velx=numpy.loadtxt("ux.dat")
    vely=numpy.loadtxt("uy.dat")
    gradx=numpy.gradient(velx)
    grady=numpy.gradient(vely)
    div=gradx[0]+grady[1]
    pylab.figure()    
    pylab.imshow(velx)
    pylab.colorbar()
    pylab.figure()
    pylab.imshow(vely)
    pylab.colorbar()
    pylab.figure()
    pylab.imshow(div)
    pylab.colorbar()
    read_poisson(velx,vely)

def calculate_untouched_fields():
    import pyximport
    pyximport.install()
    import poisson
    phase=numpy.loadtxt("phase300000.dat")
    dims=phase.shape
    
    velx=numpy.loadtxt("velocityx300000.dat")
    vely=numpy.loadtxt("velocityy300000.dat")
    
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
    #putting walls
    geometry[0,:]=0
    geometry[-1,:]=0
    
    zerolevel=numpy.where(geometry==0)
    scalar=numpy.zeros(dims)
    scalar2=numpy.zeros(dims)
    scalar[zerolevel]=0.0
    
    
    velx[negative]=0.0
    vely[negative]=0.0
    velx[zerolevel]=interface_velocity
    velx[0,:]=0.0
    vely[-1,:]=0.0
    vely[zerolevel]=0.0
    
    gradx=numpy.gradient(velx)
    grady=numpy.gradient(vely)
    div=gradx[1]+grady[0]
    div_cross=gradx[0]+grady[1]
    pylab.figure()
    pylab.imshow(div)
    pylab.figure()
    pylab.imshow(div_cross)
    
    
    domain=numpy.where(geometry==1)
    coors=numpy.array(zip(domain[0],domain[1]))
    print domain
    
    scalar=poisson.poisson_iteration(scalar,coors,div)
    pylab.figure()
    pylab.imshow(scalar)
    pylab.colorbar()
    pylab.figure()
    pylab.contour(scalar)
    pylab.colorbar()
     
    pylab.figure()
    pylab.imshow(scalar)
    pylab.colorbar()
    
    grad_scalar=numpy.gradient(scalar)
    
    fieldx=velx-grad_scalar[0]
    fieldy=vely-grad_scalar[1]
    numpy.savetxt("result_velx.dat",fieldx)
    numpy.savetxt("result_vely.dat",fieldy)
    gradx_scalar=numpy.gradient(fieldx)
    grady_scalar=numpy.gradient(fieldy)
       
    pylab.figure()
    pylab.imshow(gradx_scalar[1]+gradx_scalar[0])
    pylab.colorbar()
    pylab.figure()
    pylab.imshow(gradx[1]+grady[0])
    pylab.colorbar()
    #pylab.figure()
    #pylab.imshow(velx)
    #pylab.colorbar()
    #pylab.figure()
    #pylab.imshow(vely)
    #pylab.colorbar()
def show_divergence():
    phase=numpy.loadtxt("phase300000.dat")
    dims=phase.shape
    
    velx=numpy.loadtxt("velocityx300000.dat")
    vely=numpy.loadtxt("velocityy300000.dat")
    
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
    #putting walls
    geometry[0,:]=0
    geometry[-1,:]=0
    
    pylab.figure()
    pylab.imshow(geometry)
    pylab.title("Geometry")

    zerolevel=numpy.where(geometry==0)
    velx[negative]=interface_velocity
    vely[negative]=0.0
    velx[zerolevel]=interface_velocity
    velx[0,:]=0.0
    velx[-1,:]=0.0
    vely[zerolevel]=0.0
    
    gradx=numpy.gradient(velx)
    grady=numpy.gradient(vely)
    div=gradx[1]+grady[0]
    pylab.figure()
    pylab.imshow(div)
    pylab.colorbar()
    pylab.title("Initial divergence")
     
    velx=velx-interface_velocity
    #velx[0,:]=0.0
    #velx[-1,:]=0.0
    ux=velx.transpose()
    uy=vely.transpose()

    
    #returning back geometry
    geometry[0,:]=1
    geometry[-1,:]=1
    pylab.savetxt("ux_non.dat",ux)
    pylab.savetxt("uy_non.dat",uy)
    pylab.savetxt("geometry_non.dat",geometry.transpose(),fmt="%d")
    
    gradux=numpy.gradient(ux)
    graduy=numpy.gradient(uy)
    pylab.figure()
    div=gradux[0]+graduy[1]
    pylab.imshow(div)
    pylab.title("After shift divergence")
    pylab.colorbar()

def show_files():
    ux=numpy.loadtxt("ux_non.dat")
    uy=numpy.loadtxt("uy_non.dat")
    geom=numpy.loadtxt("geometry.dat")
   
    pylab.figure()
    pylab.imshow(ux)
    pylab.colorbar()
    
    pylab.figure()
    pylab.plot(ux[:,100])
    pylab.colorbar()
    
    pylab.figure()
    pylab.plot(ux[1000,:])
    pylab.colorbar()
    

    phase=numpy.loadtxt("phase300000.dat")
    dims=phase.shape
    
    velx=numpy.loadtxt("velocityx300000.dat")
    vely=numpy.loadtxt("velocityy300000.dat")
  
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
    print interface_velocity
    pylab.figure()
    pylab.plot(velx[:,((z1+z2)/2)%dims[1]])

def show_normals():
    ux=numpy.loadtxt("ux_non.dat")
    uy=numpy.loadtxt("uy_non.dat")
    geom=numpy.loadtxt("geometry.dat")

    zerolevel=numpy.where(geom==0)
    cx=[0,1,0,-1,0,1,-1,-1,1]
    cy=[0,0,1,0,-1,1,1,-1,-1] 
    dirs=zip(cx,cy)
    
    for x,y in zip(zerolevel[0],zerolevel[1]):
        for dirx,diry in dirs:
            posx=x+dirx
            posy=y+diry
            if geom[posx,posy]==1:
                geom[posx,posy]=2
                continue
    numpy.savetxt("geometry_adjacent.dat",geom)
    pylab.figure()
    pylab.imshow(geom)
    pylab.colorbar()
    
    #generating normals
    interface=numpy.where(geom==2)
    normx=numpy.zeros(geom.shape)
    normy=numpy.zeros(geom.shape)
    for x,y in zip(interface[0],interface[1]):
        nx=0
        ny=0
        for dirx,diry in dirs:
            posx=x+dirx
            posy=y+diry
            if geom[posx,posy]==0:
			    nx=nx+dirx
			    ny=ny+diry
        normx[posx,posy]=nx
        normy[posx,posy]=ny
	
	#bulk=numpy.where(geom==1)
	#normx[bulk]=None
	#normy[bulk]=None
	
	#void=numpy.where(geom==-1)
	#normx[void]=None
	#normy[void]=None
	
	#normx[zerolevel]=None
	#normy[zerolevel]=None	
    #M = numpy.zeros(geom.shape, dtype='bool')
    #M[bulk] = True
    #normx2 = ma.masked_array(normx, mask=M)
    #normy2 = ma.masked_array(normy, mask=M)
    #print M
	normx2=normx.transpose()
	normy2=normy.transpose()
	
    pylab.figure(figsize=(30,5))	
    Q=pylab.quiver(normx2,normy2,minlength=0.1,scale=100)
    #qk=pylab.quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',fontproperties={'weight': 'bold'})  


def show_hydro(name):
    geom=numpy.loadtxt("geometry_non.dat")
    ux_initial=numpy.loadtxt("ux_non.dat")
    uy_initial=numpy.loadtxt("uy_non.dat")
    
    density=numpy.loadtxt("density"+name)
    ux=numpy.loadtxt("velx"+name)
    uy=numpy.loadtxt("vely"+name)
    
    pylab.figure()
    pylab.imshow(geom)
    
    pylab.figure()
    pylab.imshow(ux_initial)
    pylab.colorbar()
    
    pylab.figure()
    pylab.imshow(uy_initial)
    pylab.colorbar()
    
    #pylab.figure()
    #pylab.imshow(density)
    #pylab.colorbar()
    
    pylab.figure()
    pylab.imshow(uy)
    pylab.colorbar()
    
    pylab.figure()
    pylab.title("Center")
    pylab.plot(uy[:,100])
    pylab.plot(ux_initial[:,100])
	
    pylab.figure()
    pylab.title("1000")
    pylab.plot(uy[1000,:])
    pylab.plot(ux_initial[1000,:])
    
    pylab.figure()
    pylab.title("844")
    pylab.plot(uy[844,:])
    pylab.plot(ux_initial[844,:])
    
    pylab.figure()
    pylab.imshow(ux)
    pylab.colorbar()
    
    pylab.figure()
    pylab.title("Center")
    pylab.plot(ux[:,100])
    pylab.plot(uy_initial[:,100])
    
    pylab.figure(figsize=(30,5))
    pylab.quiver(uy.transpose(),ux.transpose(),minlength=0.001)

def show_mass(name):
    mass=numpy.loadtxt("density_mass"+name)
    pylab.figure()
    pylab.imshow(mass)
    
    pylab.figure()
    pylab.plot(mass[100,:])   

def show_streamlines(name):
    from PyNGL import Ngl

    phase=numpy.loadtxt("geometry_non.dat")
    density=numpy.loadtxt("density"+name)
    ux=numpy.loadtxt("velx"+name)
    uy=numpy.loadtxt("vely"+name)
    #ux=numpy.loadtxt("uy_non.dat")
    #uy=numpy.loadtxt("ux_non.dat")
    dims=ux.shape        
 
    x,y=numpy.mgrid[0:dims[0],0:dims[1]]        

        
    positive=numpy.where(phase>0.0)
    negative=numpy.where(phase<0.0)
    
    large=numpy.where(numpy.logical_or(x<20,x>dims[0]-20))
    #print z1,z2 #,bounds
        #print y<z1-30
        #x_short=x[::3,::25]*deltay
        #y_short=y[::3,::25]*deltax
        #[numpy.where(phase_numpy>0.0)]
        #ux[negative]=None
        #ux[large]=None
        #vz_diff_mask[bounds]=None
        
    #ux_mask=ux[::100,::5]
    #uy_mask=uy[::100,::5]
    #x_mask=x[::100,::5]
    #y_mask=y[::100,::5]

        
        #pylab.figure(figsize=(20,10))
        #pylab.quiver(y_mask,x_mask,ux_mask,uy_mask,headwidth=6,minlength=0.1)
        #pylab.contour(array['phi'],[0.0],linewidths=[3])
        
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,"test")
    resources = Ngl.Resources()
     
        #uvar = file.variables["U_GRD_6_ISBL"]
        #vvar = file.variables["V_GRD_6_ISBL"]
        #if hasattr(uvar,"units"):
        #  resources.tiMainString = "GRD_6_ISBL (u,v " + uvar.units + ")"
        #else:
        #resources.tiMainString = "GRD_6_ISBL"
        #if hasattr(uvar,"_FillValue"):
        #    resources.vfMissingUValueV = uvar._FillValue
        # if hasattr(vvar,"_FillValue"):
        #    resources.vfMissingVValueV = vvar._FillValue
        
        
        #resources.tiMainFont    = "Times-Roman"
        #resources.tiMainOn=True
        #resources.tiMainString="Ca=0.22 "
           
        #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 3*0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 0.25
    #resources.wkPaperSize="A5"
    resources.nglFrame = False
    #resources.vfXArray=numpy.linspace(0.0,1.0,len(ux[1,::50]))
    #resources.vfYArray=numpy.linspace(0.0,15.0,len(ux[::5,1]))
    
       
        #resources2=Ngl.Resources()
        #resources2.tiMainFont    = "Times-Roman"
        #resources2.tiXAxisString = "streamlines"
        #resources2.tiMainOn=True
        #resources2.tiMainString="Ca=0.22"
        
        #resources2.wkPaperSize="A5"
        #resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
        #resources2.vpWidthF  = 3*0.25
        #resources2.nglFrame = False
        
        #resources2.cnLineLabelsOn = False   # Turn off contour line labels.
        #resources2.cnLinesOn      = False   # Turn off contour lines.
        #resources2.cnFillOn       = False    # Turn on contour fill.
        #resources2.cnInfoLabelOn   = False 
  
        
        #resources2.cnLevelSelectionMode = "ExplicitLevels"  # Select contour levels. 
        #resources2.cnMinLevelValF       = 0.0
        #resources2.cnMaxLevelValF       = 0.001
        
        ##resources2.cnLevelSpacingF      = 0.0
        #resources2.cnLevelCount=1
        #resources2.cnLevels=[0.0]
        ##resources2.cnLineThicknesses=[3]
        #resources2.cnMonoLineThickness=True
        #resources2.cnLineThicknessF=3.0
        
        #resources2.lbLabelBarOn=False
        #resources2.lbLabelsOn=False
        #resources2.sfXArray=numpy.linspace(0.0,15.0,len(ux[1,:]))
        #resources2.sfYArray=numpy.linspace(0.0,1.0,len(ux[:,1]))
        
        #plot = Ngl.streamline(wks,uvar[0,::2,::2],vvar[0,::2,::2],resources) 
        #print vz_diff.shape
        #print vy.shape
        #x,y=numpy.mgrid[0:dims[0],0:dims[1]]
        #vx=numpy.sin(x)*numpy.sin(y)
        #vy=numpy.cos(x)*numpy.cos(y)
    plot=Ngl.streamline(wks,ux[::5,::2],uy[::5,::2],resources)
        #Ngl.contour(wks,phase[::5,::50],resources2)        
        #Ngl.contour(wks,phase,resources2)        
        #plot=Ngl.streamline(wks,vx,vy,resources)
    Ngl.end()
        


if __name__=="__main__":
    #calculate_untouched_fields()
    #show_divergence()
    #show_files()
    #show_normals()
    #show_fields()
    name="0003900.dat"
    
    #show_hydro(name)
    #show_streamlines(name)
    show_mass(name)
    pylab.show()

