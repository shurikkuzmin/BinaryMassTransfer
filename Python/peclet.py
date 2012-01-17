#!/usr/bin/python
import numpy
import os
import subprocess
import math
import pylab

def run_simulations():
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    scale_str=[str(x) for x in [2,4,6,8,10,15,20,40]]
    for counter,dir_temp in enumerate(capillary_str):
        os.chdir(dir_temp)
        for scale in range(0,len(scale_str)-counter):
            subprocess.call(['mkdir','-p',scale_str[scale]])
           
            os.chdir(scale_str[scale])
            subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/ScalingPeclet/'+dir_temp+"/"+scale_str[scale]+"/concentration.dat","."])
            os.chdir("..")
        os.chdir("..")

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)

def analyze_particular_simulation(dir_name):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    scale_str=[str(x) for x in [2,4,6,8,10,15,20,40]]
    gamma=0.0728
    diam=1.5e-3    

    ind=capillary_str.index(dir_name)    
    os.chdir("ScalingPeclet/"+dir_name) 
    
    name_geom="geometry.dat"
    name_velx="velx0200000.dat"
    name_vely="vely0200000.dat"
    
    phase=numpy.loadtxt(name_geom)
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
    #print "Bubble length=",(z2-z1)/200.0
    #print "Slug length=",15.0-(z2-z1)/200.0
    #print "Liquid_velocity",numpy.sum(velx[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2)
    #print "Gas holdup",float(len(numpy.where(array['phi']<0)[0]))/(dims[1]*(dims[0]-2))
    #prof=phase[:,((z1+z2)/2)%dims[1]]
    #width=append(Get_Zero(prof))     
    
    interface_velocity=velx[dims[0]/2,z2%dims[1]]
    interface_velocity=vely[0,0]    
    capillary=interface_velocity*(2.0/3.0)/math.sqrt(8*0.04*0.04/9.0)
    reynolds=interface_velocity*dims[0]/(2.0/3.0)        
    phys_velocity=math.sqrt(capillary*reynolds*gamma/(1000*diam))        
    superficial_liquid=numpy.sum(velx[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2)
    gas_holdup=float(len(numpy.where(phase<0)[0]))/(dims[1]*(dims[0]-2))

    #bubble_velocities.append(phys_velocity)
     
    print "Interface velocity=",interface_velocity
    print "Capillary=",capillary
    print "Reynolds=",reynolds
    
    deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)  
        
    for scale in range(0,len(scale_str)-ind):
        os.chdir(scale_str[scale])
        conc=numpy.loadtxt("concentration.dat")
        
        new_conc=numpy.array(zip(conc[1:,0]-conc[:-1,0],conc[1:,1]-conc[:-1,1],conc[1:,2]))
        print new_conc
        arr=numpy.array(numpy.where(numpy.abs(conc[1:,2])<0.9)).ravel()
      
        ###Periodic Central analysis 
        #coefficient=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        coefficient=new_conc[arr,1]/(200*3000)/(new_conc[arr,0]*scale*(1.0-numpy.abs(new_conc[arr,2])))
        print coefficient
  
        
        #coefficient_last=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-array[1:,1]/((1-gas_holdup)*200*3000)))        
        #coefficient_aver=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-0.5*(array[1:,1]+array[0:-1,1])/((1-gas_holdup)*200*3000)))                
        #coefficient_outlet=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))        
        #coefficient_outlet_sim=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-conc_outlet))
        #coefficient_inlet_aver=inlet_new_conc[:,1]/(200*3000)/(deltat*inlet_new_conc[:,0]*(1.0-inlet_conc_aver[1:,1]/((1-gas_holdup)*200*3000)))        

        
        #coeffs_inlet.append(numpy.mean(coefficient[numpy.where(conc[1:,0]>500000)]))        
        #coeffs_aver.append(numpy.mean(coefficient_aver[numpy.where(array[1:,0]>1500000)]))
        #coeffs_outlet.append(numpy.mean(coefficient_outlet[numpy.where(array[1:,0]>3500000)]))
        #coeffs_open.append(numpy.mean(coefficient_open))
                
        pylab.figure(1)
        #pylab.plot(conc[1:,0],coefficient) #,styles[counter]) 
        pylab.plot(coefficient)
        pylab.title(r'''$Ca='''+str(capillary)[0:4]+'''$''')
        #pylab.xlim(xmax=600000)
        #pylab.ylim(ymax=0.00001)
        
        pylab.figure(2)
        pylab.plot(conc[:,0],conc[:,2])
        #pylab.plot(open_conc[:,0],coefficient_open,styles[counter])        
        #pylab.plot(open_time,coefficient_open,styles[counter])        

        #pylab.figure(6)
        #pylab.plot(inlet_conc_aver[1:,0],coefficient_inlet_aver,styles[counter])        
        
        os.chdir("..")
    pylab.figure(1)    
    pylab.savefig("scalingpeclet"+str(capillary)[0]+str(capillary)[2:4]+".eps",format="EPS",dpi=300)

def check_particular_contours(dir_name,factor,subtract):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    capillaries=[0.097,0.254,0.526,0.750,1.040]
    scale=[2,4,8,10,20,40]
    scale_str=[str(x) for x in [2,4,8,10,20,40]]
    os.chdir("ScalingPeclet/"+dir_name)

    c_levels=numpy.arange(0.2,1.0,0.2)
    counter=capillary_str.index(dir_name) 
    capillary=capillaries[counter]
    for scale_counter in range(0,len(scale_str)-subtract):
        #subprocess.call(['mkdir','-p',scale_str[scale]])
        os.chdir(scale_str[scale_counter])
        name=str(20000*factor*scale[len(scale_str)-subtract-1]/scale[scale_counter])
        name="density"+"0"*(7-len(name))+name+".dat"        
        print "Filename=",name
        print "Scale=",scale_str[scale_counter]        
        #subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/ScalingPeclet/'+dir_temp+"/"+scale_str[scale]+"/concentration.dat","."])
        arr=numpy.loadtxt(name)
        dims=arr.shape
        gas_holdup=float(len(numpy.where(arr>0.0)[0]))/(dims[1]*(dims[0]-2))

        print "Overall concentration=",numpy.sum(arr[numpy.where(arr>0.0)])/(gas_holdup*dims[1]*(dims[0]-2))       
        pylab.figure(scale_counter)
        pylab.imshow(arr)
        pylab.colorbar()
        
        
        pylab.figure(99,figsize=(12,3))
        
        cont=pylab.contour(arr,levels=c_levels,linestyles="dotted",colors=['k'],linewidths=[2],extent=(0,15,0,1))
        if scale_counter==len(scale_str)-subtract-1:        
            pylab.clabel(cont)
        #else:
        #    cont.clabel(fmt=scale_str[scale_counter],fontsize="12",inline_spacing=10)
        os.chdir("..")
        
    fig=pylab.figure(99)
    
    fig.subplots_adjust(bottom=0.2,top=0.8)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$x$''',fontsize=30)
    pylab.ylabel(r'''$y$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    pylab.savefig("contourlines_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)
    
if __name__=="__main__":
    #modify_file()
    #analyze_particular_simulation("42")
    #check_particular_contours("9",1,0)    
    #check_particular_contours("21",4,1)
    #check_particular_contours("42",5,4)
    check_particular_contours("60",10,5)    
    #check_particular_contours("84",10,5)    
    pylab.show()
