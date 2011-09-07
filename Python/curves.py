#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy
import pylab
import math
import os
import subprocess

def Produce_Bunch():
    print os.getcwd()
    #capillary_str=["3","5","8","10","20","40"] #,"60","80"]
    #bubble_velocities=numpy.array([0.023,0.043,0.073,0.059,0.202,0.437]) #,0.671,0.902])
    
    capillary_str=["3","5","8","10","20","40","60","80"]
    bubble_velocities=numpy.array([0.023,0.043,0.073,0.059,0.202,0.437,0.671,0.902])

    styles=['kv','ks','ko','k^','k>','k<','kD','kp']
    aver_coefficient=[]
    fig=pylab.figure(99)
    for counter,dir_temp in enumerate(capillary_str):
        dir_name="Results/"+dir_temp
        name=dir_name+"/concentration.dat"
        array=numpy.loadtxt(name)
        #print array
        new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        print "Data=",new_conc
        coefficient=new_conc[:,1]/(200*3000)/(4.64e-7*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        pylab.figure(99)
        pylab.plot(array[1:,0],coefficient,styles[counter])
        aver_coefficient.append(numpy.mean(coefficient[len(coefficient)/2:]))
    legs=[r'''$U_{\mathrm{bubble}}='''+str(vel)[:4]+r'''$''' for vel in bubble_velocities]
    pylab.xlabel('Iterations',fontsize=20)
    pylab.ylabel(r'''$k_L a, \mathrm{s^{-1}}$''',fontsize=20)
    pylab.legend(legs,fancybox=True,labelspacing=0.1)
    pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    print "Aver_coefficient=",aver_coefficient
    return aver_coefficient
def Produce_Pictures():
    print os.getcwd()
    capillary_str=["3","5","8","10","20","40","60","80"]
    capillaries=numpy.array([0.026,0.047,0.080,0.065,0.222,0.479,0.736,0.989])
    
    #capillary=capillaries[capillary_str.index(dir_name)]
    for counter,dir_name in enumerate(capillary_str):
        dir_name="Results/"+dir_name
        os.chdir(dir_name)

        capillary=capillaries[counter]
        file_list=[]    
        for root,dirs,files in os.walk(os.getcwd()):
            for file in files:
                if file[0:7]=="density":
                    file_list.append(file)
        #for file_name in sorted(file_list):
        #    arr=numpy.loadtxt(file_name)
        #    print file_name        
        #    pylab.imshow(arr,origin="lower")
        #    pylab.title(r'''$Ca='''+str(capillary)+r'''$''')
        #    pylab.text(1500, 100, "Time="+file_name[7:-4])
        #    pylab.savefig("map"+file_name[7:-4]+".png",format="PNG")
        #    pylab.clf()
        subprocess.call(["mencoder","mf://*.png","-mf","fps=5:type=png","-ovc","lavc","-lavcopts","vcodec=mpeg4","-o","output.mp4"])
        #subprocess.call(["mencoder","mf://*.png","-mf","fps=5:type=png","-ovc","lavc","-lavcopts","vcodec=wmv1","-o","ca0"+str(capillary)[2:]+".avi"])
        os.chdir("../..")
    #   #print file
    #    read_file=numpy.loadtxt(file)
    #    flux.append(read_file[0].tolist())
    
def Produce_Curves():
    capillaries=numpy.array([0.026,0.047,0.080,0.065,0.222,0.479,0.736,0.989])
    bubble_lengths_non=numpy.array([5.225,5.22,5.215,4.965,5.25,5.565,5.51,5.505])
    slug_lengths_non=numpy.array([9.775,9.78,9.785,10.035,9.75,9.435,9.49,9.495])
    bubble_velocities=numpy.array([0.023,0.043,0.073,0.059,0.202,0.437,0.671,0.902])
    holdups=numpy.array([0.3064,0.2935,0.2800,0.2667,0.2532,0.2494,0.2361,0.2300])
    superficial_liq_lbm=numpy.array([0.00135,0.00238,0.00383,0.00317,0.00970,0.01959,0.02903,0.03841])
    widths=numpy.array([0.04015,0.05819,0.08552,0.07622,0.12295,0.15123,0.16457,0.17233])
    velocities_lbm=numpy.array([0.00148,0.00270,0.00454,0.00371,0.01256,0.02713,0.04164,0.05597]) 
    reynolds=numpy.array([0.449486,0.82096,1.37800,1.12643,3.80721,8.22257,12.61778,16.96062])
    a=numpy.array([2.55389217042,2.54312761748,2.53266472281,2.41422758192,2.52531357555,2.65187070664,2.62692244618,2.63404014475])
    diameter=1.5e-3
    hydraulic=2*diameter
    #hydraulic=diameter    
    
    bubble_lengths=diameter*bubble_lengths_non
    slug_lengths=diameter*slug_lengths_non
    bubble_diameters=diameter*(1.0-2.0*widths)

    superficial_gas=bubble_velocities*holdups    
    superficial_liquid=superficial_liq_lbm*bubble_velocities/velocities_lbm
    
    diffusion=1.118e-9
    deltat=4.64e-7
    print bubble_velocities/velocities_lbm    
    print "Liquid superficial=",superficial_liquid
    print "Bubble_velocities=",bubble_velocities
    print "Gas_superficial=",superficial_gas
    print "Iterations=",diameter*diameter*widths*widths/(8.0*diffusion*deltat)    
    fourier=diffusion*(bubble_lengths-diameter)/(bubble_velocities*widths*widths*diameter*diameter)
    aver_coeff=Produce_Bunch()

    yue=2.0/hydraulic*numpy.sqrt(diffusion*bubble_velocities/(bubble_lengths+slug_lengths))*numpy.power(bubble_lengths/(bubble_lengths+slug_lengths),0.3)
    bercic=0.111*math.sqrt(diffusion/2e-9)*numpy.power(superficial_gas+superficial_liquid,1.19)/numpy.power((1-holdups)*(bubble_lengths+slug_lengths),0.57)    
    vanbaten=2.0/math.sqrt(math.pi)*numpy.sqrt(diffusion*bubble_velocities/(bubble_lengths-hydraulic))*4\
            *(bubble_lengths-hydraulic)/(hydraulic*(bubble_lengths+slug_lengths))\
            +2.0*math.sqrt(2.0)/math.pi*numpy.sqrt(diffusion*bubble_velocities/hydraulic)*4/(bubble_lengths+slug_lengths)
    vanbaten_my=4.0*numpy.sqrt(diffusion*bubble_velocities/math.pi)*numpy.sqrt(bubble_lengths-bubble_diameters)/((bubble_lengths+slug_lengths)*diameter)\
               +2.0*math.sqrt(2.0)*numpy.sqrt(diffusion*bubble_velocities)*numpy.sqrt(bubble_diameters)/((bubble_lengths+slug_lengths)*diameter)
    
    fig=pylab.figure()    
    pylab.plot(bubble_velocities,yue,'k^',markersize=8)
    pylab.plot(bubble_velocities,bercic,'ko',markersize=8)
    pylab.plot(bubble_velocities,vanbaten,'ks',markersize=8)
    pylab.plot(bubble_velocities,aver_coeff,'kd',markersize=8)
    pylab.plot(bubble_velocities,vanbaten_my,'kh',markersize=8)
    leg=["Yue","Bercic","Van Baten","Simulations","Adjusted Van Baten"]
    pylab.legend(leg,loc=2)
    pylab.xlabel(r'''$U_{\mathrm{bubble}},\,\mathrm{m/s}$''',fontsize=30)
    pylab.ylabel(r'''$k_L a,\,\mathrm{s^{-1}}$''',fontsize=30)
    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    fig.subplots_adjust(left=0.17,bottom=0.17) 
    pylab.savefig("comparison_correlations.eps",format="EPS",dpi=300)
    print "Fourier=",fourier
    
    pylab.show() 
def Produce_Mass_Bunch_Sailfish():
    #physical parameters
    gamma=0.0728
    diam=1.5e-3    
    print os.getcwd()
 
    styles=['kv','ks','ko','k^','k>','k<','kD','kh']
    aver_coefficient=[]
    bubble_velocities=[]
    dirs=[3,5,8,10,20,40,60,80]
    #dirs=[3,5,8,10,20]    
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
        
        bubble_velocities.append(phys_velocity)
        
        print "Interface velocity=",interface_velocity
        print "Capillary=",capillary
        print "Reynolds=",reynolds
  
        dir_name="SailfishMassBGK/"+str(dir_name)+"/"
        #dir_name="SailfishMassJos/"+str(dir_name)+"/"     
        name=dir_name+"concentration.dat"
        array=numpy.loadtxt(name)
        new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)
            
        print "Data=",new_conc
        coefficient=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        pylab.figure(99)
        pylab.plot(array[1:,0],coefficient,styles[counter])
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
    pylab.legend(legs,fancybox=True,labelspacing=0.1)
    #pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    pylab.figure()
    print phys_velocity
    pylab.plot(bubble_velocities,aver_coefficient,'ro')
    print "Aver_coefficient=",aver_coefficient
    return aver_coefficient


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
        dir_name="SailfishMassJos/"+str(dir_name)+"/"     
        name=dir_name+"concentration.dat"
        array=numpy.loadtxt(name)
        #new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        #deltat=interface_velocity/phys_velocity*diam/(dims[0]-2)   
        #print "Data=",new_conc
        #coefficient=new_conc[:,1]/(200*3000)/(deltat*new_conc[:,0]*(1.0-numpy.abs(new_conc[:,2])))
        coefficient=(superficial_liquid+gas_holdup*interface_velocity)*phys_velocity/interface_velocity\
                   /(15*diam)*numpy.log((1-numpy.abs(array[:,2]))/(1-numpy.abs(array[:,3])))
        pylab.figure(99)
        pylab.plot(array[:,0],coefficient,styles[counter])
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
    pylab.legend(legs,fancybox=True,labelspacing=0.1)
    #pylab.savefig("steady_state.eps",format="EPS",dpi=300)
    
    pylab.figure()
    print phys_velocity
    pylab.plot(bubble_velocities,aver_coefficient,'ro')
    print "Aver_coefficient=",aver_coefficient
    return aver_coefficient



if __name__=="__main__":
    #Produce_Pictures()
    #Produce_Curves()
    Produce_Mass_Bunch_Sailfish()
    #Produce_Mass_Bunch_Jos()    
    pylab.show()