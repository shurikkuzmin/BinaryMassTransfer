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
    dirs=[3,5,8,10,20,40,60,80]
    #dirs=[3,5,8,10]
    coeffs_inlet_aver=[]
    coeffs_aver=[]
    coeffs_outlet=[]
    coeffs_open=[]
   
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
      
        deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)
      
        ###Periodic Inlet analysis
  
        inlet_dir_name="SailfishMassBGK/"+str(dir_name)+"/"
        inlet_name=inlet_dir_name+"concentration.dat"
        inlet_conc_aver=numpy.loadtxt(inlet_name)
        inlet_new_conc=numpy.array(zip(inlet_conc_aver[1:,0]-inlet_conc_aver[:-1,0],inlet_conc_aver[1:,1]-inlet_conc_aver[:-1,1],inlet_conc_aver[1:,2]))
        print inlet_conc_aver.shape
        print inlet_new_conc.shape
        
      
        ###Periodic Central analysis 
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
        coefficient_inlet_aver=inlet_new_conc[:,1]/(200*3000)/(deltat*inlet_new_conc[:,0]*(1.0-inlet_conc_aver[1:,1]/((1-gas_holdup)*200*3000)))        

        
        coeffs_inlet_aver.append(numpy.mean(coefficient_inlet_aver[numpy.where(inlet_conc_aver[1:,0]>1000000)]))        
        coeffs_aver.append(numpy.mean(coefficient_aver[numpy.where(array[1:,0]>1500000)]))
        coeffs_outlet.append(numpy.mean(coefficient_outlet[numpy.where(array[1:,0]>3500000)]))
        coeffs_open.append(numpy.mean(coefficient_open))
                
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

        pylab.figure(6)
        pylab.plot(inlet_conc_aver[1:,0],coefficient_inlet_aver,styles[counter])        

        #pylab.title("Open boundaries")
    #pylab.legend(legs,fancybox=True,labelspacing=0.1)
    #pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    print "Coeffs_aver=",coeffs_aver
    print "Coeffs_outlet=",coeffs_outlet
    print "Coeffs_open=",coeffs_open
    print "Coeff_inlet_aver=",coeffs_inlet_aver
    
    #print "Aver_coefficient=",aver_coefficient
    pylab.figure(1)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)
    pylab.ylim(ymax=3.0)
    pylab.savefig("steady_state_per_last.eps",dpi=300)
    
    pylab.figure(2)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.savefig("steady_state_per_aver.eps",dpi=300)
    pylab.ylim(ymax=3.0)

    pylab.figure(3)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.ylim(ymax=3.0)
    pylab.savefig("steady_state_per_middle_outlet.eps",dpi=300)

    
    pylab.figure(4)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.ylim(ymax=3.0)

    pylab.figure(5)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities[0:4]]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.xlim(xmin=0)    
    pylab.savefig("steady_state_jos_inlet_outlet.eps",dpi=300)
    
    pylab.figure(6)
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)   
    pylab.ylim(ymax=3.0)
    pylab.savefig("steady_state_inlet_aver.eps",dpi=300)
    
    #pylab.ylim(ymax=3.0,ymin=-2.0)
   #pylab.savefig("steady_state_average.eps",format="EPS",dpi=300)
    
    #pylab.figure()
    #print phys_velocity
    
    #pylab.plot(bubble_velocities,aver_coefficient,'ro')
    #print "Aver_coefficient=",aver_coefficient
    #return aver_coefficient

def Produce_Double():

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
 
 
    coeffs_overall=[]
    coeffs_first=[]
    coeffs_second=[]
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
  
        results_dir_name="SailfishMassCentralDouble/"+str(dir_name)+"/"     
        
        vx=numpy.loadtxt("SailfishResultsCenter/"+str(dir_name)+"/ux.dat").transpose()
        coefficient_overall=[]
        coefficient_first=[]
        coefficient_second=[]
        for file_number in numpy.arange(1160000,4000000,280000):
            density=numpy.loadtxt(results_dir_name+"density_jos_par"+str(file_number)+".dat")
            inlet_conc=numpy.sum(density[1:-1,1]*vx[1:-1,1])/numpy.sum(vx[1:-1,1])
            outlet_conc=numpy.sum(density[1:-1,-2]*vx[1:-1,-2])/numpy.sum(vx[1:-1,-2])
            middle_conc=numpy.sum(density[1:-1,density.shape[1]/2]*vx[1:-1,-2])/numpy.sum(vx[1:-1,-2])
            coefficient_overall.append((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
                   /(2*15*diam)*numpy.log((1-numpy.abs(inlet_conc))/(1-numpy.abs(outlet_conc))))
            coefficient_first.append((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
                   /(15*diam)*numpy.log((1-numpy.abs(inlet_conc))/(1-numpy.abs(middle_conc))))
            coefficient_second.append((superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
                   /(15*diam)*numpy.log((1-numpy.abs(middle_conc))/(1-numpy.abs(outlet_conc))))
            
        coefficient_overall=numpy.array(coefficient_overall)
        coefficient_first=numpy.array(coefficient_first)
        coefficient_second=numpy.array(coefficient_second)
        time_iterations=numpy.arange(1160000,4000000,280000)
        inds=numpy.where(time_iterations>3300000)        
        print coefficient_overall[inds[0]]        
        coeffs_overall.append(numpy.mean(coefficient_overall[inds[0]]))
        coeffs_first.append(numpy.mean(coefficient_first[inds[0]]))
        coeffs_second.append(numpy.mean(coefficient_second[inds[0]]))
        #numpy.savetxt(results_dir_name+"coefficient.dat",zip(time_iterations,coefficient))           
        
        pylab.figure(1)
        pylab.plot(time_iterations,coefficient_overall,styles[counter])        
        pylab.figure(2)
        pylab.plot(time_iterations,coefficient_first,styles[counter])        
        pylab.figure(3)
        pylab.plot(time_iterations,coefficient_second,styles[counter])        
        
        #aver_coefficient.append(numpy.mean(coefficient[len(coefficient)/2:]))
    print "Coeffs_overall=", coeffs_overall
    print "Coeffs_first=",coeffs_first
    print "Coeffs_second=",coeffs_second
    file_names=["steady_double_overall.eps","steady_double_first.eps","steady_double_second.eps"]        
    for i in range(1,4):    
        pylab.figure(i)        
        legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
        pylab.xlabel('Iterations',fontsize=20)
        pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
        pylab.legend(legs,fancybox=True,labelspacing=0.1)
        pylab.ylim(ymax=9)
        pylab.savefig(file_names[i-1],dpi=300)
       
    #pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    #print "Aver_coefficient=",aver_coefficient



if __name__=="__main__":
    Produce_Double()    
    pylab.show()