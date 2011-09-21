#!/usr/bin/python
import numpy
import pylab

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
    
    domain=numpy.where(geometry==1)
    coors=zip(domain[0],domain[1])
    print domain
    
    for counter in range(1,100):
        for i,j in coors:
            topj=(j+1+dims[1])%dims[1]
            bottomj=(j-1+dims[1])%dims[1]
            delta=0.25*(scalar[i,topj]+scalar[i,bottomj]+scalar[i+1,j]+scalar[i-1,j]-4*scalar[i,j]+4*div[i,j])
            scalar2[i,j]=scalar[i,j]+delta
        scalar=scalar2
        print counter

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
    pylab.figure()
    pylab.imshow(velx)
    pylab.colorbar()
    pylab.figure()
    pylab.imshow(vely)
    pylab.colorbar()

if __name__=="__main__":
    calculate_untouched_fields()
    pylab.show()

