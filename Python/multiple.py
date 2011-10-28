#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy
import pylab
import math
import os
import subprocess
import glob


def Produce_Single():
    #physical parameters
    gamma=0.0728
    diam=1.5e-3    
    print os.getcwd()
 
    styles=['kv','ks','ko','k^','k>','k<','kD','kh']
    aver_coefficient=[]
    bubble_velocities=[]
    #dirs=[3,5,8,10,20,40,60,80]
    #dirs=[3,5,8,10,20,40,60,80]
    dirs=[3,5,8,10]    
    
   
    for counter,dir_name in enumerate(dirs):    
        data_dir_name="SailfishData/"+str(dir_name)+"/"
        file_name=data_dir_name+"capillary200000.npz"        
        arr=numpy.load(file_name)        
        phase=arr['phi']
        velx=arr['v'][0]
        vely=arr['v'][1]
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
        superficial_liquid=numpy.sum(velx[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2)
        gas_holdup=float(len(numpy.where(phase<0)[0]))/(dims[1]*(dims[0]-2))

        bubble_velocities.append(phys_velocity)
        
        print "Interface velocity=",interface_velocity
        print "Capillary=",capillary
        print "Reynolds=",reynolds
  
        ###Periodic analysis 
        per_dir_name="SailfishMassCentralPeriodic/"+str(dir_name)+"/"
        name=per_dir_name+"concentration.dat"
        density_name=per_dir_name+"density3960000.dat"
        ux_name=per_dir_name+"ux.dat"
        
        array=numpy.loadtxt(name)
        density=numpy.loadtxt(density_name)
        #density=density.transpose()
        ux=numpy.loadtxt(ux_name)
        ux=ux.transpose()
        print "Density=",density.shape
        print "ux=",ux.shape
        conc_outlet=numpy.sum(ux[1:-1,-1]*density[1:-1,-1])/numpy.sum(ux[1:-1,-1])
        print "Outlet=",conc_outlet        
        new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)
        
        ###Open boundaries options
        open_dir_name="SailfishMassCentralSingle/"+str(dir_name)+"/"     
        open_conc=numpy.loadtxt(open_dir_name+"concentration.dat")       
         
        coefficient_open=[]        
        for file_number in numpy.arange(1160000,4000000,280000):
            density_number=numpy.loadtxt(open_dir_name+"density"+str(file_number)+".dat")
            inlet_conc=numpy.sum(density_number[1:-1,1]*ux[1:-1,1])/numpy.sum(ux[1:-1,1])
            outlet_conc=numpy.sum(density_number[1:-1,-2]*ux[1:-1,-2])/numpy.sum(ux[1:-1,-2])
            coefficient_open.append((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
                   /(15*diam)*(1-abs(inlet_conc))/(1-abs(outlet_conc)))
        open_time=numpy.arange(1160000,4000000,280000) 
        print coefficient_open        
        #files_names=glob.glob(results_dir_name+"density*.dat")
        #vx=numpy.loadtxt("SailfishResults/"+str(dir_name)+"/ux.dat").transpose()
        #files_names=sorted(files_names)
        #for file_example in files_names:
        #    print file_example,
        #    density=numpy.loadtxt(file_example)
        #    inlet_conc=numpy.sum(density[1:-1,1]*vx[1:-1,1])/numpy.sum(vx[1:-1,1])
        #    outlet_conc=numpy.sum(density[1:-1,-2]*vx[1:-1,-2])/numpy.sum(vx[1:-1,-2])
        #    coefficient.append((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
        #           /(15*diam)*numpy.log((1-numpy.abs(inlet_conc))/(1-numpy.abs(outlet_conc))))
        #    time_iterations.append(int(file_example[-11:-4]))
        #numpy.savetxt(results_dir_name+"coefficient.dat",zip(time_iterations,coefficient))

        #coefficient_open=((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
        #           /(15*diam)*numpy.log((1-numpy.abs(open_conc[:,2]))/(1-numpy.abs(open_conc[:,3]))))        

        
        #coefficient=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        coefficient_last=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-array[1:,1]/((1-gas_holdup)*200*3000)))        
        coefficient_aver=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-0.5*(array[1:,1]+array[0:-1,1])/((1-gas_holdup)*200*3000)))                
        coefficient_outlet=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))        
        coefficient_outlet_sim=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-conc_outlet))
        pylab.figure(1)
        pylab.plot(array[1:,0],coefficient_last,styles[counter])        
        pylab.figure(2)
        pylab.plot(array[1:,0],coefficient_aver,styles[counter])
        pylab.figure(3)
        pylab.plot(array[1:,0],coefficient_outlet,styles[counter])
        pylab.figure(4)
        pylab.plot(array[1:,0],coefficient_outlet_sim,styles[counter])
        pylab.figure(5)
        #pylab.plot(open_conc[:,0],coefficient_open,styles[counter])        
        pylab.plot(open_time,coefficient_open,styles[counter])        
        #pylab.title("Open boundaries")
    #pylab.legend(legs,fancybox=True,labelspacing=0.1)
    #pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    #print "Aver_coefficient=",aver_coefficient
    pylab.figure(1)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)
    pylab.ylim(ymax=3.0)
    pylab.figure(2)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.ylim(ymax=3.0)

    pylab.figure(3)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.ylim(ymax=3.0)

    pylab.figure(4)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.ylim(ymax=3.0)

    pylab.figure(5)
    pylab.xlim(xmin=0)    
    #pylab.ylim(ymax=3.0,ymin=-2.0)
   #pylab.savefig("steady_state_average.eps",format="EPS",dpi=300)
    
    #pylab.figure()
    #print phys_velocity
    
    #pylab.plot(bubble_velocities,aver_coefficient,'ro')
    #print "Aver_coefficient=",aver_coefficient
    #return aver_coefficient
    
def Produce_Mass_Bunch_Jos():

    #physical parameters
    gamma=0.0728
    diam=1.5e-3    
    print os.getcwd()
 
    styles=['kv','ks','ko','k^','k>','k<']
    aver_coefficient=[]
    bubble_velocities=[]
    interface_velocities=[]
    liquid_velocities=[]
    gas_velocities=[]
    holdups=[]    
    #dirs=[3,5,8,10,20,40]
    dirs=[3,5,8,10]    
    fig=pylab.figure(99)

    for counter,dir_name in enumerate(dirs):    
        data_dir_name="SailfishData/"+str(dir_name)+"/"
        file_name=data_dir_name+"capillary200000.npz"        
        arr=numpy.load(file_name)        
        phase=arr['phi']
        velx=arr['v'][0]
        vely=arr['v'][1]
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
        superficial_liquid=numpy.sum(velx[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2)
        gas_holdup=float(len(numpy.where(phase<0)[0]))/(dims[1]*(dims[0]-2))
        
        interface_velocities.append(interface_velocity)
        bubble_velocities.append(phys_velocity)
        liquid_velocities.append(superficial_liquid)
        gas_velocities.append(gas_holdup*interface_velocity)        
                
        holdups.append(gas_holdup)
        
        print "Interface velocity=",interface_velocity
        print "Capillary=",capillary
        print "Reynolds=",reynolds
  
        #dir_name="SailfishMassBGK/"+str(dir_name)+"/"
        
        results_dir_name="SailfishMassJos/"+str(dir_name)+"/"     
        
        files_names=glob.glob(results_dir_name+"density*.dat")
        vx=numpy.loadtxt("SailfishResults/"+str(dir_name)+"/ux.dat").transpose()
        files_names=sorted(files_names)
        coefficient=[]
        time_iterations=[]
        for file_example in files_names:
            print file_example,
            density=numpy.loadtxt(file_example)
            inlet_conc=numpy.sum(density[1:-1,1]*vx[1:-1,1])/numpy.sum(vx[1:-1,1])
            outlet_conc=numpy.sum(density[1:-1,-2]*vx[1:-1,-2])/numpy.sum(vx[1:-1,-2])
            coefficient.append((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
                   /(15*diam)*numpy.log((1-numpy.abs(inlet_conc))/(1-numpy.abs(outlet_conc))))
            time_iterations.append(int(file_example[-11:-4]))
        numpy.savetxt(results_dir_name+"coefficient.dat",zip(time_iterations,coefficient))

            
        #array=numpy.loadtxt(results_dir_name+"concentration.dat")        
        #new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        #deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)   
        #print "Data=",new_conc
        #coefficient=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        
        #Calculation of coefficients
        #inlet_conc=numpy.sum(density[1:-1,1]*vx[1:-1,1])/numpy.sum(vx[1:-1,1])
        #outlet_conc=numpy.sum(density[1:-1,-2]*vx[1:-1,-2])/numpy.sum(vx[1:-1,-2])
        #print "Inlet conc=",inlet_conc
        #print "Outlet conc=",outlet_conc
        
        #coefficient=(superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
        #           /(15*diam)*numpy.log((1-numpy.abs(inlet_conc))/(1-numpy.abs(outlet_conc)))

        pylab.figure(99)
        #pylab.plot(array[:,0],coefficient,styles[counter])
        pylab.plot(time_iterations,coefficient,styles[counter])        
        aver_coefficient.append(numpy.mean(coefficient[len(coefficient)/2:]))
        
        
        #coefficient=(superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
        #           /(15*diam)*numpy.log((1-numpy.abs(array[:,2]))/(1-numpy.abs(array[:,3])))
        #pylab.figure(99)
        #pylab.plot(array[:,0],coefficient,styles[counter])
        #aver_coefficient.append(numpy.mean(coefficient[len(coefficient)/2:]))
    
            
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    #pylab.legend(legs,fancybox=True,labelspacing=0.1)
    pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    print "Aver_coefficient=",aver_coefficient
    #pylab.xlim(0,1000000)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)
    #pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    pylab.figure()
    print phys_velocity
    pylab.plot(bubble_velocities,aver_coefficient,'ro')
    print "Aver_coefficient=",aver_coefficient
    return aver_coefficient

def Produce_Mass_Bunch_Jos_Grid():

    #physical parameters
    gamma=0.0728
    diam=1.5e-3    
    print os.getcwd()
 
    styles=['kv','ks','ko','k^','k>','k<']
    aver_coefficient=[]
    bubble_velocities=[]
    interface_velocities=[]
    liquid_velocities=[]
    gas_velocities=[]
    holdups=[]    
    #dirs=[3,5,8,10,20,40]
    dirs=[5,8,10]    
    fig=pylab.figure(99)

    for counter,dir_name in enumerate(dirs):    
        data_dir_name="SailfishData/"+str(dir_name)+"/"
        file_name=data_dir_name+"capillary200000.npz"        
        arr=numpy.load(file_name)        
        phase=arr['phi']
        velx=arr['v'][0]
        vely=arr['v'][1]
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
        superficial_liquid=numpy.sum(velx[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2)
        gas_holdup=float(len(numpy.where(phase<0)[0]))/(dims[1]*(dims[0]-2))
        
        interface_velocities.append(interface_velocity)
        bubble_velocities.append(phys_velocity)
        liquid_velocities.append(superficial_liquid)
        gas_velocities.append(gas_holdup*interface_velocity)        
                
        holdups.append(gas_holdup)
        
        print "Interface velocity=",interface_velocity
        print "Capillary=",capillary
        print "Reynolds=",reynolds
  
        #dir_name="SailfishMassBGK/"+str(dir_name)+"/"
        results_dir_name="SailfishMassJosDouble/"+str(dir_name)+"/"     
        
        files_names=glob.glob(results_dir_name+"density*.dat")
        vx=numpy.loadtxt("SailfishResults/"+str(dir_name)+"/ux.dat").transpose()
        files_names=sorted(files_names)
        coefficient=[]
        time_iterations=[]
        print files_names
        for file_example in files_names:
            print file_example,
            density=numpy.loadtxt(file_example)
            inlet_conc=numpy.sum(density[1:-1,1]*vx[1:-1,1])/numpy.sum(vx[1:-1,1])
            outlet_conc=numpy.sum(density[1:-1,-2]*vx[1:-1,-2])/numpy.sum(vx[1:-1,-2])
            coefficient.append((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
                   /(30*diam)*numpy.log((1-numpy.abs(inlet_conc))/(1-numpy.abs(outlet_conc))))
            time_iterations.append(int(file_example[-11:-4]))

        numpy.savetxt(results_dir_name+"coefficient.dat",zip(time_iterations,coefficient))           
        
        #array=numpy.loadtxt(results_dir_name+"concentration.dat")        
        #new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        #deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)   
        #print "Data=",new_conc
        #coefficient=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        
        #Calculation of coefficients
        #inlet_conc=numpy.sum(density[1:-1,1]*vx[1:-1,1])/numpy.sum(vx[1:-1,1])
        #outlet_conc=numpy.sum(density[1:-1,-2]*vx[1:-1,-2])/numpy.sum(vx[1:-1,-2])
        #print "Inlet conc=",inlet_conc
        #print "Outlet conc=",outlet_conc
        
        #coefficient=(superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
        #           /(15*diam)*numpy.log((1-numpy.abs(inlet_conc))/(1-numpy.abs(outlet_conc)))

        pylab.figure(99)
        #pylab.plot(array[:,0],coefficient,styles[counter])
        pylab.plot(time_iterations,coefficient,styles[counter])        
        aver_coefficient.append(numpy.mean(coefficient[len(coefficient)/2:]))
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    #pylab.legend(legs,fancybox=True,labelspacing=0.1)
    pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    print "Aver_coefficient=",aver_coefficient
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1,loc=4)
    pylab.savefig("steady_state_double_jos.eps",format="EPS",dpi=300)
    
    #pylab.figure()
    #print phys_velocity
    #pylab.plot(bubble_velocities,aver_coefficient,'ro')
    #print "Aver_coefficient=",aver_coefficient
    #return aver_coefficient




if __name__=="__main__":
    Produce_Single()    
    pylab.show()