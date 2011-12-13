#!/usr/bin/python
import numpy
import pylab
import scipy.special
import math

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

def film_analytical():
    num_terms=30
    c0=0.0
    cs=1.0
    u_bubble=0.05
    y,x=0.01*numpy.mgrid[1:101,1:501]
    
    c=numpy.zeros_like(x)
    diffusion=1.0/3.0*(1-0.5)
    
    sign=-1
    for i in range(0,num_terms):
        sign=-sign        
        c=c+sign*(scipy.special.erfc((y+2*i)/numpy.sqrt(4*x*diffusion/u_bubble))+scipy.special.erfc((2*(i+1)-y)/numpy.sqrt(4*x*diffusion/u_bubble)))
    print "X=",x
    print "Y=",y
    print "C=",c
    c=c0-(c0-cs)*c
    pylab.imshow(c)
    pylab.colorbar()

def compare_film(file_dir):
    concentration=numpy.loadtxt(file_dir+"film0050000.dat")
    uy=numpy.loadtxt(file_dir+"uy_initial.dat")
    
    dims=concentration.shape    
    print dims
    pylab.figure()    
    pylab.imshow(concentration[:,:-4],extent=(0,10,0,1))
    pylab.title("Simulations")
    pylab.colorbar()    
    #pylab.figure()
    #pylab.plot(concentration[:,dims[1]/2])
    #pylab.figure()
    #pylab.plot(concentration[:,0])

    num_terms=30
    c0=0.5
    cs=1.0
    u_bubble=0.05*20   
    y,x=0.01*numpy.mgrid[1:101,1:1001]
    
    c=numpy.zeros_like(x)
    diffusion=1.0/3.0*(1.0/1.4-0.5)
    
    res=y/numpy.sqrt(4*x*diffusion/u_bubble)
    print res.shape
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        c=c+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x*diffusion/u_bubble))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x*diffusion/u_bubble)))
        
    #print "X=",x
    #print "Y=",y
    #print "C=",c
    c=c0-(c0-cs)*c
    pylab.figure()    
    pylab.imshow(c,extent=(0,10,0,1))
    pylab.colorbar()
    pylab.title("Analytics")

    pylab.figure()
    pylab.imshow(uy)
    pylab.colorbar()
    
    pylab.figure()
    pylab.plot(concentration[0,:])
    pylab.figure()
    pylab.plot(concentration[:,0])
    pylab.figure()
    pylab.plot(concentration[:,dims[1]-1])
    
    c_levels=numpy.arange(0.0,1.0,0.1)
    print c_levels
    pylab.figure(figsize=(10,1))
    c1=pylab.contour(concentration[:,:-4],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(c,levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)    
    
    print numpy.where(concentration>1.8)

def compare_three_films(file_dir):
    conc_antibb =numpy.loadtxt(file_dir+"film0050000.dat")
    conc_inamuro=numpy.loadtxt(file_dir+"film_inamuro0050000.dat")
    conc_outflow=numpy.loadtxt(file_dir+"film_outflow0050000.dat")    
    
    pylab.figure()    
    pylab.imshow(conc_antibb[:,:-20],extent=(0,10,0,1))
    pylab.title("AntiBB simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_inamuro[:,:-20],extent=(0,10,0,1))
    pylab.title("Inamuro simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_outflow[:,:-20],extent=(0,10,0,1))
    pylab.title("Partial Inamuro simulations")
    pylab.colorbar()

    c_levels=numpy.arange(0.6,1.0,0.1)
    print c_levels
    pylab.figure(99,figsize=(10,1))
    c1=pylab.contour(conc_antibb[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(conc_inamuro[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)
    c3=pylab.contour(conc_outflow[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))



    num_terms=30
    c0=0.5
    cs=1.0
    u_bubble=0.05*20   
    y,x=0.01*numpy.mgrid[1:101,1:1001]
    
    c=numpy.zeros_like(x)
    diffusion=1.0/3.0*(1.0/1.4-0.5)
    
    res=y/numpy.sqrt(4*x*diffusion/u_bubble)
    print res.shape
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        c=c+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x*diffusion/u_bubble))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x*diffusion/u_bubble)))
        
    c=c0-(c0-cs)*c
    pylab.figure()    
    pylab.imshow(c,extent=(0,10,0,1))
    pylab.colorbar()
    pylab.title("Analytics")

    pylab.figure(99,)
    c_anal=pylab.contour(c[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c_anal,fontsize=9, inline=1)    


def compare_full_profiles(file_dir):
    conc_antibb =numpy.loadtxt(file_dir+"FullProfile/"+"film_antibb0050000.dat")
    conc_inamuro=numpy.loadtxt(file_dir+"FullProfile/"+"film_outflow0050000.dat")
    
    print "Shape=",conc_antibb.shape
    pylab.figure()    
    pylab.imshow(conc_antibb,extent=(0,20,0,1))
    pylab.title("AntiBB simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_inamuro,extent=(0,20,0,1))
    pylab.title("Inamuro simulations")
    pylab.colorbar()    

    c_levels=numpy.arange(0.2,1.0,0.1)
    print c_levels
    pylab.figure(99,figsize=(20,1))
    c1=pylab.contour(conc_antibb,levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(conc_inamuro,levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)



    num_terms=100
    c0=0.0
    cs=1.0
    ny=40
    omega=1.8
    u_bubble=0.05
    diffusion=1.0/3.0*(1.0/omega-0.5)
    pe=u_bubble*ny/diffusion
    print "Pe=",pe
    y,x=0.005*numpy.mgrid[1:201,1:8001]
    #pe=pe/2
    c=numpy.zeros_like(x)
    
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        c=c+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x/pe))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x/pe)))
    c=c0-(c0-cs)*c
    #c=scipy.special.erfc(y/numpy.sqrt(4.0*x/pe))    
    print c    
    pylab.figure()    
    pylab.imshow(c,extent=(0,20,0,0.5))
    pylab.colorbar()
    pylab.title("Analytics")

    pylab.figure(99)
    #pylab.figure(figsize=(10,1))
    c_anal=pylab.contour(c,levels=c_levels,linestyles="dashed",extent=(0.0,10.0,0.0,0.5))
    #pylab.clabel(c_anal,fontsize=9, inline=1)

def film_reconstruction(file_dir):
    #bessel_zeros=scipy.special.jv(-3.0/4.0,num_terms)
    wmfivetwo=numpy.array([0.209965, 0.0181809, 0.0069193, 0.00373185, 0.00236737, 0.00165045, \
                           0.00122407, 0.000948389, 0.000759097, 0.000623069, 0.000521773, \
                           0.000444147, 0.000383242, 0.000334505, 0.000294845, 0.000262104, \
                           0.000234734, 0.000211602, 0.000191861, 0.000174867, 0.000160124, \
                           0.000147245, 0.000135921, 0.000125909, 0.000117008, 0.000109058, \
                           0.000101925, 0.0000954983, 0.0000896863, 0.0000844115, 0.0000796084, \
                           0.0000752211, 0.0000712022, 0.0000675107, 0.0000641113, 0.0000609733, \
                        0.0000580702, 0.0000553784, 0.0000528777, 0.00005055])
    wmcubic=numpy.array([0.137974, 0.0368546, 0.0213316, 0.0150107, 0.0115796, 0.00942524, \
                         0.00794677, 0.00686924, 0.00604903, 0.0054038, 0.00488295, \
                         0.00445368, 0.00409379, 0.00378771, 0.00352422, 0.003295, 0.00309378, \
                         0.00291572, 0.00275704, 0.00261474, 0.00248641, 0.00237009, \
                         0.00226416, 0.0021673, 0.00207838, 0.00199648, 0.00192078, \
                         0.00185062, 0.0017854, 0.00172462, 0.00166784, 0.00161468, \
                         0.00156481, 0.00151792, 0.00147376, 0.0014321, 0.00139273, \
                         0.00135547, 0.00132015, 0.00128662])
    
    bessel_zeros=numpy.array([1.454997085, 2.927133004, 3.857578101, 4.601777732, 5.240824067, \
                  5.809787864, 6.327694301, 6.806248002, 7.253261690, 7.674259300, \
                  8.073317953, 8.453548974, 8.817390869, 9.166797024, 9.503360989, \
                  9.828402986, 10.14303138, 10.44818742, 10.74467855, 11.03320360, \
                  11.31437223, 11.58872006, 11.85672070, 12.11879537, 12.37532065, \
                  12.62663484, 12.87304320, 13.11482231, 13.35222370, 13.58547689, \
                  13.81479204, 14.04036213, 14.26236488, 14.48096438, 14.69631251, \
                  14.90855018, 15.11780841, 15.32420926, 15.52786670, 15.72888729, \
                  15.92737088, 16.12341118, 16.31709625, 16.50850900, 16.69772758, \
                  16.88482577, 17.06987328, 17.25293611, 17.43407678, 17.61335459, \
                  17.79082587, 17.96654416, 18.14056039, 18.31292309, 18.48367852, \
                  18.65287083, 18.82054216, 18.98673282, 19.15148136, 19.31482468, \
                  19.47679813, 19.63743561, 19.79676966, 19.95483148, 20.11165108, \
                  20.26725729, 20.42167785, 20.57493946, 20.72706782, 20.87808772, \
                  21.02802303, 21.17689679, 21.32473123, 21.47154783, 21.61736732, \
                  21.76220974, 21.90609448, 22.04904029, 22.19106530, 22.33218710, \
                  22.47242269, 22.61178856, 22.75030070, 22.88797461, 23.02482533, \
                  23.16086744, 23.29611511, 23.43058210, 23.56428178, 23.69722713, \
                  23.82943078, 23.96090501, 24.09166175, 24.22171263, 24.35106896, \
                  24.47974174, 24.60774171, 24.73507930, 24.86176469, 24.98780781])
    y,x=0.01*numpy.mgrid[1:101,1:2001]
    c=numpy.zeros_like(x)
    cerf=numpy.zeros_like(x)
    
    num_terms=40
    cwall=1
    c0=0
    coeffm=wmfivetwo/wmcubic*(c0-cwall)
    omega=1.8
    u_bubble=0.1*20
    diffusion=1.0/3.0*(1.0/omega-0.5)

    pe=u_bubble/diffusion
    
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        cerf=cerf+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x/pe))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x/pe)))
    cerf=c0-(c0-cwall)*cerf
    
    for i in range(0,num_terms):
        c=c+coeffm[i]*scipy.special.jv(0.25,bessel_zeros[i]*bessel_zeros[i]*0.5*y*y)*numpy.sqrt(y)*numpy.exp(-(bessel_zeros[i]**4)*x*diffusion/u_bubble)
    c=c+cwall    
    pylab.figure()
    pylab.imshow(c,extent=(0,10,0,0.5))
    pylab.colorbar()
 
    c_levels=numpy.arange(0.2,1.0,0.1)
 
    conc_antibb =numpy.loadtxt(file_dir+"FullProfile/"+"film_antibb0050000.dat")
    conc_inamuro=numpy.loadtxt(file_dir+"FullProfile/"+"film_outflow0050000.dat")
    dims=conc_antibb.shape    
    print "Dims=",dims
    print "Pe=",pe
    
    pylab.figure()    
    pylab.imshow(conc_antibb[0:dims[0]/2,:],extent=(0,20,0,0.5))
    pylab.title("AntiBB simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_inamuro[0:dims[0]/2,:],extent=(0,20,0,0.5))
    pylab.title("Inamuro simulations")
    pylab.colorbar()
    
    pylab.figure(99,figsize=(20,1))
    c1=pylab.contour(conc_antibb,levels=c_levels,extent=(0.0,20.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(conc_inamuro,levels=c_levels,extent=(0.0,20.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)
    
    c_anal=pylab.contour(c,levels=c_levels,linestyles="dashed",extent=(0.0,10.0,0.0,0.5))
    #c_anal_erf=pylab.contour(cerf,levels=c_levels,linestyles="dotted",extent=(0.0,10.0,0.0,0.5))


if __name__=="__main__":
    #file_name="../Benchmarks/density0001000.dat"
    file_name="../Benchmarks/conc_initial.dat"    
    file_dir="../Benchmarks/"
    #read_file(file_name)
    
    #read_files(file_dir)
    #show_bessel()    
    #read_film(file_dir)    
    #film_analytical()    
    #compare_film(file_dir)
    #compare_three_films(file_dir)
    #compare_full_profiles(file_dir)    
    film_reconstruction(file_dir)
    pylab.show()
