#!/usr/bin/python
import numpy
import os
import subprocess
import math
import pylab
import scipy
import scipy.optimize

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

def check_average_concentration(dir_name,factor,subtract):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    capillaries=[0.097,0.254,0.526,0.750,1.040]
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    velocities_slug=[0.0046,0.0108,0.0209,0.0293,0.0397]
    velocities_gas=[0.0016,0.0041,0.0080,0.0101,0.0135]
    


    scale=[2,4,8,10,20,40]
    #scale=[8,10,20,40]
    scale_str=[str(x) for x in scale]
    os.chdir("ScalingPeclet/"+dir_name)
    
    geometry=numpy.loadtxt("geometry.dat")
    gas_holdup=len(numpy.where(geometry>0)[0])/(3000.0*200.0)
    counter=capillary_str.index(dir_name) 
    capillary=capillaries[counter]
   
    #styles=["k+","ko","ks","k^","kv","kd"]    
    
    styles=['k-','k--','k-.','k:','k.','k,']    
    
    exam_points=[]
    legs=[]
    for scale_counter in range(0,len(scale_str)-subtract):
        os.chdir(scale_str[scale_counter])
        #exam_point=int(20*factor*scale[len(scale_str)-subtract-1]/scale[scale_counter])
        #exam_point=int(24.0/(velocities[counter]*scale[scale_counter]))        
        exam_point=10;        
        #name=str(20000*factor*scale[len(scale_str)-subtract-1]/scale[scale_counter])
        #name="density"+"0"*(7-len(name))+name+".dat"        
        #print "Filename=",name
        #print "Scale=",scale_str[scale_counter]        
        #subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/ScalingPeclet/'+dir_temp+"/"+scale_str[scale]+"/concentration.dat","."])
        #arr=numpy.loadtxt(name)
        #dims=arr.shape
        #gas_holdup=float(len(numpy.where(arr>0.0)[0]))/(dims[1]*(dims[0]-2))

        #print "Overall concentration=",numpy.sum(arr[numpy.where(arr>0.0)])/(gas_holdup*dims[1]*(dims[0]-2))       
       # pylab.figure(scale_counter)
        #pylab.imshow(arr)
        #pylab.colorbar()
        concentration=numpy.loadtxt("concentration.dat")
        aver_concentration=concentration[1:,1]/(3000*200*gas_holdup)
        pylab.figure(1)
        pylab.plot(velocities[counter]*scale[scale_counter]*concentration[1:,0],aver_concentration,styles[scale_counter])
        pylab.figure(2)
        
        #this plot is in terms of unit cells        
        #pylab.plot(velocities[counter]*scale[scale_counter]*concentration[1:,0]/3000.0,3000.0/(velocities[counter]*scale[scale_counter]*concentration[1:,0])*numpy.log(1.0/(1.0-aver_concentration)),styles[scale_counter],linewidth=2)                
        #this plot is in terms of time        
        pylab.plot(concentration[1:,0]/1000.0,3000.0/(velocities[counter]*scale[scale_counter]*concentration[1:,0])*numpy.log(1.0/(1.0-aver_concentration)),styles[scale_counter],linewidth=2)                

        pylab.figure(99)
        pylab.plot(velocities[counter]*scale[scale_counter]*concentration[1:,0]/3000.0,3000.0/(velocities[counter]*scale[scale_counter]*concentration[1:,0])*numpy.log(1.0/(1.0-concentration[1:,2])),styles[scale_counter],linewidth=2)                
      

        pylab.figure(3)
        window_len=10     
        conc=[]
        for i in range(0,len(aver_concentration)-window_len):
            conc.append(3000.0/(velocities[counter]*scale[scale_counter]*1000.0*window_len)\
                        *numpy.log((1.0-aver_concentration[i])/(1.0-aver_concentration[i+window_len])))   
        pylab.plot(conc,styles[scale_counter])
        pylab.ylim(ymin=0.0,ymax=0.5)
        pylab.figure(4)
        vanbaten=[]
        mult_vanbaten=(velocities_gas[counter]+velocities_slug[counter])*scale[scale_counter]
        for i in range(0,len(aver_concentration)-window_len):
            vanbaten.append((aver_concentration[i+window_len]-aver_concentration[i])*3000.0/(1000*mult_vanbaten*window_len*(1.0-aver_concentration[i+window_len/2])))   
        pylab.plot(velocities[counter]*scale[scale_counter]*concentration[0:len(aver_concentration)-window_len,0]/3000.0,\
                   vanbaten,styles[scale_counter])
        pylab.ylim(ymin=0.0,ymax=0.5)
        
        pylab.figure(100)
        vanbaten_final=[]
        mult_vanbaten=(velocities_gas[counter]+velocities_slug[counter])*scale[scale_counter]
        for i in range(0,len(aver_concentration)-window_len):
            vanbaten_final.append((aver_concentration[i+window_len]-aver_concentration[i])*3000.0/(1000*mult_vanbaten*window_len*(1.0-concentration[i+window_len/2,2])))   
        pylab.plot(velocities[counter]*scale[scale_counter]*concentration[0:len(aver_concentration)-window_len,0]/3000.0,vanbaten_final,styles[scale_counter])
        pylab.ylim(ymin=0.0,ymax=1.0)
        pylab.xlim(xmax=15)
        
        
        #pylab.figure(99,figsize=(12,3))
        print 3000.0/(velocities[counter]*scale[scale_counter]*concentration[exam_point,0])*numpy.log(1.0/(1.0-aver_concentration[exam_point-1]))
        legs.append(scale_str[scale_counter])        
        #cont=pylab.contour(arr,levels=c_levels,linestyles="dotted",colors=['k'],linewidths=[2],extent=(0,15,0,1))
        #if scale_counter==len(scale_str)-subtract-1:        
        #    pylab.clabel(cont)
        #else:
        #    cont.clabel(fmt=scale_str[scale_counter],fontsize="12",inline_spacing=10)
        os.chdir("..")
        
    #fig=pylab.figure(99)
    
    fig=pylab.figure(2)    
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.18)
    pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.5,ymin=0.0)
    #pylab.xlabel(r'''$\mathrm{Cell\,Units}$''',fontsize=30)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    
    pylab.ylabel(r'''$k_L a \frac{L_{\mathrm{unit}}}{U_{\mathrm{liq}}+U_{\mathrm{gas}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    #pylab.savefig("aver_conc_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)
    pylab.savefig("aver_conc_scale_ca_time"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)


    fig=pylab.figure(99)    
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.18)
    pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.75,ymin=0.0)
    pylab.xlabel(r'''$\mathrm{Cell\,Units}$''',fontsize=30)
    #pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    
    pylab.ylabel(r'''$k_L a \frac{L_{\mathrm{unit}}}{U_{\mathrm{liq}}+U_{\mathrm{gas}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    #pylab.savefig("aver_conc_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)
    pylab.savefig("aver_vanbaten_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)


    fig=pylab.figure(4)    
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.18)
    pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.75,ymin=0.0)
    pylab.xlim(xmax=15)
    pylab.xlabel(r'''$\mathrm{Cell\,Units}$''',fontsize=30)
    #pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    
    pylab.ylabel(r'''$k_L a \frac{L_{\mathrm{unit}}}{U_{\mathrm{liq}}+U_{\mathrm{gas}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    #pylab.savefig("aver_conc_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)
    pylab.savefig("vanbaten_aver_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)

    fig=pylab.figure(100)    
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.18)
    pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.75,ymin=0.0)
    pylab.xlabel(r'''$\mathrm{Cell\,Units}$''',fontsize=30)
    #pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    
    pylab.ylabel(r'''$k_L a \frac{L_{\mathrm{unit}}}{U_{\mathrm{liq}}+U_{\mathrm{gas}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    #pylab.savefig("aver_conc_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)
    pylab.savefig("vanbaten_full_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)




def check_average_jos_concentration(dir_name,scales):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    capillaries=[0.097,0.254,0.526,0.750,1.040]
    
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    velocities_slug=[0.0046,0.0108,0.0209,0.0293,0.0397]
    velocities_gas=[0.0016,0.0041,0.0080,0.0101,0.0135]
    
    os.chdir("JosScalingPeclet/"+dir_name)
    
    geometry=numpy.loadtxt("geometry.dat")
    gas_holdup=len(numpy.where(geometry>0)[0])/(3000.0*200.0)
    counter=capillary_str.index(dir_name) 
    capillary=capillaries[counter]
   
    #styles=["k+","ko","ks","k^","kv","kd"]    

    styles=['k-','k--','k-.','k:','k.','k,','ko','ks']    
    
    exam_points=[]
    legs=[]
    for scale_counter,scale in enumerate(scales):
        os.chdir(scale)
        #arr=numpy.loadtxt(name)
        #dims=arr.shape
        #gas_holdup=float(len(numpy.where(arr>0.0)[0]))/(dims[1]*(dims[0]-2))
        #print "Overall concentration=",numpy.sum(arr[numpy.where(arr>0.0)])/(gas_holdup*dims[1]*(dims[0]-2))       
        #pylab.figure(scale_counter)
        #pylab.imshow(arr)
        #pylab.colorbar()
        concentration=numpy.loadtxt("concentration.dat")
        dims=concentration.shape
 
       
        
        #exam_point=int(20*factor*scale[len(scale_str)-subtract-1]/scale[scale_counter])
        #exam_point=int(24.0/(velocities[counter]*float(scale)))        
        exam_point=dims[0]/2        
        aver_concentration=concentration[:,1]/(3000*200*gas_holdup)
        
        window_len=10     
        conc=[]
        for i in range(0,dims[0]-window_len):
            conc.append(3000.0/(velocities[counter]*float(scale)*1000*window_len)\
                        *numpy.log((1.0-aver_concentration[i])/(1.0-aver_concentration[i+window_len])))   
           
        
        pylab.figure(1)
        pylab.plot(velocities[counter]*float(scale)*concentration[:,0],aver_concentration,styles[scale_counter])
        pylab.figure(2)
        pylab.plot(velocities[counter]*float(scale)*concentration[:,0]/3000.0,3000.0/(velocities[counter]*float(scale)*concentration[:,0])*numpy.log(1.0/(1.0-aver_concentration)),styles[scale_counter],linewidth=2)                
        #pylab.figure(99,figsize=(12,3))
        
        pylab.figure(3)
        pylab.plot(conc,styles[scale_counter])
        print 3000.0/(velocities[counter]*float(scale)*concentration[exam_point,0])*numpy.log(1.0/(1.0-aver_concentration[exam_point-1]))
        legs.append(scale)
        
        pylab.figure(4)
        conc_inlet_outlet=numpy.log((1.0-concentration[:,2])/(1.0-concentration[:,3])) #*(float(scale)*velocities_slug[counter])*1000        
        #pylab.plot(conc_inlet_outlet)
        #pylab.plot((concentration[1:,1]-concentration[:-1,1]))
        multiplier=float(scale)*velocities_slug[counter]/(velocities_slug[counter]+velocities_gas[counter])
        mult_vanbaten=float(scale)*(velocities_slug[counter]+velocities_gas[counter])        
        vanbaten=[]        
        for i in range(0,dims[0]-window_len):
            vanbaten.append((aver_concentration[i+window_len]-aver_concentration[i])*3000.0/(1000*mult_vanbaten*window_len*(1.0-aver_concentration[i+window_len/2])))   
        
        #pylab.plot((concentration[:,3]-concentration[:,2])*multiplier/(1.0-aver_concentration)) #/3000.0)        
        pylab.plot(vanbaten)
        pylab.plot()        
        os.chdir("..")
        
    fig=pylab.figure(3)
    pylab.ylim(ymin=0.0,ymax=1.0)
    pylab.legend(legs)   
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.05,ymin=0.0)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    pylab.ylabel(r'''$k_L a \frac{L}{U_{\mathrm{gas}}+U_{\mathrm{bubble}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.2)
    pylab.savefig("jos_aver_moving_window_ca"+str(capillary)[0:1]+str(capillary)[2:]+".eps",format="EPS",dpi=300)
    
    fig=pylab.figure(2)    
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.2)
    pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.05,ymin=0.0)
    pylab.xlabel(r'''$L$''',fontsize=30)
    pylab.ylabel(r'''$k_L a \frac{L}{U_{\mathrm{gas}}+U_{\mathrm{bubble}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    pylab.savefig("jos_aver_conc_scale_ca"+str(capillary)[0:1]+str(capillary)[2:]+".eps",format="EPS",dpi=300)

def check_average_sym_concentration(dir_name,scales):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    capillaries=[0.097,0.254,0.526,0.750,1.040]
    
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]

    os.chdir("JosScalingPeclet/"+dir_name)
    
    geometry=numpy.loadtxt("geometry.dat")
    gas_holdup=len(numpy.where(geometry>0)[0])/(3000.0*200.0)
    counter=capillary_str.index(dir_name) 
    capillary=capillaries[counter]
   
    #styles=["k+","ko","ks","k^","kv","kd"]    

    styles=['k-','k--','k-.','k:','k.','k,','ko','ks']    
    
    exam_points=[]
    legs=[]
    for scale_counter,scale in enumerate(scales):
        os.chdir(scale+"Symmetry")
        #arr=numpy.loadtxt(name)
        #dims=arr.shape
        #gas_holdup=float(len(numpy.where(arr>0.0)[0]))/(dims[1]*(dims[0]-2))
        #print "Overall concentration=",numpy.sum(arr[numpy.where(arr>0.0)])/(gas_holdup*dims[1]*(dims[0]-2))       
        #pylab.figure(scale_counter)
        #pylab.imshow(arr)
        #pylab.colorbar()
        concentration=numpy.loadtxt("concentration.dat")
        dims=concentration.shape
 
       
        
        #exam_point=int(20*factor*scale[len(scale_str)-subtract-1]/scale[scale_counter])
        #exam_point=int(24.0/(velocities[counter]*float(scale)))        
        exam_point=dims[0]/2        
        aver_concentration=concentration[:,1]/(3000*200*gas_holdup)
        
        window_len=10     
        conc=[]
        for i in range(0,dims[0]-window_len):
            conc.append(3000.0/(velocities[counter]*float(scale)*1000*window_len)\
                        *numpy.log((1.0-aver_concentration[i])/(1.0-aver_concentration[i+window_len])))   
           
        
        pylab.figure(1)
        pylab.plot(velocities[counter]*float(scale)*concentration[:,0],aver_concentration,styles[scale_counter])
        pylab.figure(2)
        pylab.plot(velocities[counter]*float(scale)*concentration[:,0]/3000.0,3000.0/(velocities[counter]*float(scale)*concentration[:,0])*numpy.log(1.0/(1.0-aver_concentration)),styles[scale_counter],linewidth=2)                
        #pylab.figure(99,figsize=(12,3))
        
        pylab.figure(3)
        pylab.plot(conc,styles[scale_counter])
        print 3000.0/(velocities[counter]*float(scale)*concentration[exam_point,0])*numpy.log(1.0/(1.0-aver_concentration[exam_point-1]))
        legs.append(scale)        
        os.chdir("..")
        
    fig=pylab.figure(3)
    pylab.ylim(ymin=0.0,ymax=1.0)
    pylab.legend(legs)   
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.05,ymin=0.0)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    pylab.ylabel(r'''$k_L a \frac{L}{U_{\mathrm{gas}}+U_{\mathrm{bubble}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.2)
    pylab.savefig("sym_aver_moving_window_ca"+str(capillary)[0:1]+str(capillary)[2:]+".eps",format="EPS",dpi=300)
    
    fig=pylab.figure(2)    
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.2)
    pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.05,ymin=0.0)
    pylab.xlabel(r'''$L$''',fontsize=30)
    pylab.ylabel(r'''$k_L a \frac{L}{U_{\mathrm{gas}}+U_{\mathrm{bubble}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    pylab.savefig("sym_aver_conc_scale_ca"+str(capillary)[0:1]+str(capillary)[2:]+".eps",format="EPS",dpi=300)


def give_overall_characteristics():
    capillary_str=["9","21","42","60","84"]
    
    capillaries_orig=[]
    reynolds_orig=[]
    velocities_orig=[]
    holdups_orig=[]
    bubble_lengths_orig=[]
    slug_lengths_orig=[]
    widths_orig=[]
    gas_velocities_orig=[]
    liq_velocities_orig=[]
    for dir_name in capillary_str:
        file_dir="HydroResultsRight/"+dir_name+"/"        
        phase=numpy.loadtxt(file_dir+"phase300000.dat")
        velx=numpy.loadtxt(file_dir+"velocityx300000.dat")
        vely=numpy.loadtxt(file_dir+"velocityy300000.dat")
        dims=phase.shape
        center=phase[dims[0]/2,:]
   
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+dims[1]
            z1=numpy.max(numpy.where(center>0.0))
        
        prof=phase[:,((z1+z2)/2)%dims[1]]
        widths_orig.append(Get_Zero(prof))     
        bubble_lengths_orig.append((z2-z1)/float(dims[0]-2))
        slug_lengths_orig.append(15.0-(z2-z1)/float(dims[0]-2))
        velocities_orig.append(velx[dims[0]/2,z2%dims[1]])
        holdups_orig.append(float(len(numpy.where(phase<0)[0]))/(dims[1]*(dims[0]-2)))        
        liq_velocities_orig.append(numpy.sum(velx[1:-1,((z1+z2+dims[1])/2)%dims[1]])/(dims[0]-2))
        gas_velocities_orig.append(holdups_orig[-1]*velocities_orig[-1])
        capillaries_orig.append(velocities_orig[-1]*(2.0/3.0)/math.sqrt(8*0.04*0.04/9.0))
        reynolds_orig.append(velocities_orig[-1]*(dims[0]-2)/(2.0/3.0))
    
    print "Bubble velocities=",velocities_orig
    print "Liquid velocities=",liq_velocities_orig
    print "Gas velocities=",gas_velocities_orig
    print "Holdups=",holdups_orig
    print "Widths=",widths_orig
    print "Capillaries=",capillaries_orig
    print "Reynolds=",reynolds_orig
    print "Bubble lengths=",bubble_lengths_orig
    print "Slug lengths=",slug_lengths_orig
    
    capillaries=[]
    reynolds=[]
    velocities=[]
    holdups=[]
    bubble_lengths=[]
    slug_lengths=[]
    widths=[]
    gas_velocities=[]
    liq_velocities=[]
 
    for counter,dir_name in enumerate(capillary_str):
        os.chdir("UnitCells/"+dir_name)
        geometry=numpy.loadtxt("geometry.dat").transpose()
        ux=numpy.loadtxt("vely0200000.dat").transpose()
        uy=numpy.loadtxt("velx0200000.dat").transpose()
        dims=geometry.shape
       
        center=geometry[dims[0]/2,:]
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+dims[1]
            z1=numpy.max(numpy.where(center>0.0))
       
        holdups.append(len(numpy.where(geometry<0)[0])/float((dims[0]-2)*dims[1]))
        velocities.append(ux[0,0])
        liq_velocities.append(ux[0,0]-numpy.sum(ux[1:-1,((z2+z1+dims[1])/2)%dims[1]])/float(dims[0]-2))
        gas_velocities.append(velocities[-1]*holdups[-1])

        prof=geometry[:,((z1+z2)/2)%dims[1]]
        widths.append(Get_Zero(prof))     
        bubble_lengths.append((z2-z1)/float(dims[0]-2))
        slug_lengths.append(15.0-(z2-z1)/float(dims[0]-2))
        capillaries.append(velocities[-1]*(2.0/3.0)/math.sqrt(8*0.04*0.04/9.0))
        reynolds.append(velocities[-1]*(dims[0]-2)/(2.0/3.0))

     
        os.chdir("../..")
    
    print "Bubble velocities=",velocities
    print "Liquid velocities=",liq_velocities
    print "Gas velocities=",gas_velocities
    print "Holdups=",holdups
    print "Widths=",widths
    print "Capillaries=",capillaries
    print "Reynolds=",reynolds_orig
    print "Bubble lengths=",bubble_lengths
    print "Slug lengths=",slug_lengths

    


def check_dependance():
    aver_coefficients=numpy.array([0.21,0.14,0.095,0.074,0.0601])
    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    bubble_lengths=numpy.array([1159,1224,1238,1193,1118])
    omega=1.99
    diffusion=1.0/3.0*(1.0/omega-0.5)
    
    
    #aver_coefficients=aver_coefficients*numpy.power(bubble_lengths/3000.0,0.3)
#    capillary_str=[str(x) for x in [9,21,42,60,84]]
#    for dir_name in capillary_str:
#        os.chdir("UnitCells/"+dir_name)
#        geometry=numpy.loadtxt("geometry.dat").transpose()
#        dims=geometry.shape
#        center=geometry[dims[0]/2,:]
#   
#        z1 = numpy.min(numpy.where(center < 0.0))
#        z2 = numpy.max(numpy.where(center < 0.0))
#        if z1==0:
#            z2=numpy.min(numpy.where(center>0.0))+dims[1]
#            z1=numpy.max(numpy.where(center>0.0))
#        print z1,z2
#        print "Bubble length=",z2-z1       
#        os.chdir("../..")

    peclets=velocities*200.0/diffusion
    peclets_range=numpy.arange(peclets[0],14000,100)
    #approximation=aver_coefficients[0]*numpy.sqrt(peclets[0])/numpy.sqrt(peclets_range)
    pylab.plot(peclets,aver_coefficients,"ko",markersize=6)
    #pylab.plot(peclets_range,approximation)
    
    #fitting procedure
    fitfunc = lambda p, x: p[0]*(x**p[1]) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    p0 = [aver_coefficients[0]*numpy.sqrt(peclets[0]), 0.5] # Initial guess for the parameters
    p1, success = scipy.optimize.leastsq(errfunc, p0[:],args=(peclets,aver_coefficients))
    
    print p1
    pylab.plot(peclets_range,p1[0]*numpy.power(peclets_range,p1[1]),"k--")
    pylab.legend(["Simulations","Correlation"])
    #pylab.text(2000,10,r'''$'''+str(p1[0])+r'''Pe^'''+str(p1[1])+r'''$''',bbox=True)
    pylab.savefig("volumetric_mass_peclet.eps",format="EPS",dpi=300)
    
    
def check_yue_dependance():
    aver_coefficients=numpy.array([0.21,0.14,0.095,0.074,0.0601])
    holdups=numpy.array([0.302466666667,0.288423333333,0.271146666667,0.252866666667,0.229446666667])

    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    velocities_liq=numpy.array([ 0.00459266431535,0.0107317293541,0.0208964048893,0.0292353942998,0.0362306511905])    
    velocities_gas=holdups*velocities    
    bubble_lengths=numpy.array([1159,1224,1238,1193,1118])
  
    omega=1.99
    diffusion=1.0/3.0*(1.0/omega-0.5)
  
   
    #aver_coefficients=aver_coefficients*numpy.power(bubble_lengths/3000.0,0.3)
    #capillary_str=[str(x) for x in [9,21,42,60,84]]
    #for counter,dir_name in enumerate(capillary_str):
    #    os.chdir("UnitCells/"+dir_name)
    #    geometry=numpy.loadtxt("geometry.dat").transpose()
    #    ux=numpy.loadtxt("vely0200000.dat").transpose()
    #    uy=numpy.loadtxt("velx0200000.dat").transpose()
    #    holdup=len(numpy.where(geometry<0)[0])/(200.0*3000.0)
    #    dims=geometry.shape
    #    center=geometry[dims[0]/2,:]
   
    #    z1 = numpy.min(numpy.where(center < 0.0))
    #    z2 = numpy.max(numpy.where(center < 0.0))
    #    if z1==0:
    #        z2=numpy.min(numpy.where(center>0.0))+dims[1]
    #        z1=numpy.max(numpy.where(center>0.0))
    #    print z1,z2
    #    print "Bubble length=",z2-z1
    #    u_liq=velocities[counter]-numpy.sum(ux[1:-1,((z2+z1+dims[1])/2)%dims[1]])/200
    #    u_gas=holdup*velocities[counter]        
    #    print "U_liq=",u_liq
    #    print "U_bubble=",velocities[counter]        
    #   print "U_gas=",u_gas
    #    print "Whole=",u_gas+u_liq
    #    print "Holdup=",holdup        
    #    os.chdir("../..")


    peclets=velocities*200.0/diffusion
    peclets_range=numpy.arange(peclets[0],14000,100)
    
    hydro_diam=200
    yue_correlation=2.0*3000.0/hydro_diam*numpy.power((diffusion/(3000.0*(velocities_gas+velocities_liq))),0.5)\
        *numpy.power((velocities/(velocities_gas+velocities_liq)),0.5)\
        *numpy.power(bubble_lengths/3000.0,0.3)
    yue2_correlation=2.0*3000.0/hydro_diam*numpy.power(1.0/(15.*peclets),0.5)\
        *numpy.power((velocities/(velocities_gas+velocities_liq)),0.5)\
        *numpy.power(bubble_lengths/3000.0,0.3)
    #fitting procedure
    fitfunc = lambda p, x: p[0]*(x**p[1]) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    p0 = [aver_coefficients[0]*numpy.sqrt(peclets[0]), -0.5] # Initial guess for the parameters
    p1, success = scipy.optimize.leastsq(errfunc, p0[:],args=(peclets,aver_coefficients))
    peclets_range=numpy.arange(peclets[0],14000,100)
    print p1
    
    fig=pylab.figure()
    pylab.plot(peclets_range,p1[0]*numpy.power(peclets_range,p1[1]),"k--")


    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    pylab.plot(peclets,aver_coefficients,"ko-",markersize=6)
    pylab.plot(peclets,yue_correlation,"ks-",markersize=6)
    pylab.xlabel(r'''$\mathrm{Peclet}$''',fontsize=30)
    pylab.ylabel(r'''$k_L a\frac{L}{U_{\mathrm{gas}}+U_{\mathrm{liq}}}$''',fontsize=30)
    pylab.legend(["Fitting","Simulations","Yue correlation"],fancybox=True)
    fig.subplots_adjust(left=0.2,bottom=0.15)    
    pylab.savefig("correlations_comparison.eps",format="EPS",dpi=300)
    

        #pylab.legend(["Simulations","Correlation"])
    #pylab.text(2000,10,r'''$'''+str(p1[0])+r'''Pe^'''+str(p1[1])+r'''$''',bbox=True)
    #pylab.savefig("volumetric_mass_peclet.eps",format="EPS",dpi=300)    

def check_vanbaten_dependance():
    aver_coefficients=numpy.array([0.21,0.14,0.095,0.074,0.0601])
    holdups=numpy.array([0.302466666667,0.288423333333,0.271146666667,0.252866666667,0.229446666667])
    widths=numpy.array([0.092,0.132,0.157,0.167,0.177])
    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    velocities_liq=numpy.array([ 0.00459266431535,0.0107317293541,0.0208964048893,0.0292353942998,0.0362306511905])    
    velocities_gas=holdups*velocities    
    bubble_lengths=numpy.array([1159,1224,1238,1193,1118])
    omega=1.99
    diffusion=1.0/3.0*(1.0/omega-0.5)


    peclets=velocities*200.0/diffusion
    
    vanbaten=3000.0/(velocities_liq+velocities_gas)*4.0*numpy.sqrt(diffusion*velocities/numpy.pi)\
             *numpy.sqrt(bubble_lengths-200.0*(1.0-2.0*widths))/(3000.0*200.0)\
             +3000.0/(velocities_liq+velocities_gas)*2.0*numpy.sqrt(2.0)*numpy.sqrt(diffusion*velocities)\
             *numpy.sqrt(200.0*(1.0-2.0*widths))/(3000.0*200.0)
    pylab.plot(peclets,aver_coefficients,"+-")
    pylab.plot(peclets,vanbaten)
    pylab.ylim(ymin=0.0)
    
if __name__=="__main__":
    #modify_file()
    #analyze_particular_simulation("42")
    #check_particular_contours("9",1,0)    
    #check_particular_contours("21",4,1)
    #check_particular_contours("42",5,4)
    #check_particular_contours("60",10,5)    
    #check_particular_contours("84",10,5)    
    
    #check_average_concentration("9",1.9,0)
    #check_average_concentration("21",4,1)
    #check_average_concentration("42",5,4)    
    #check_average_concentration("60",10,5)
    #check_average_concentration("84",10,5)    
    
    #check_average_jos_concentration("9",["2","4","6","8","10","15","20","40"])
    #check_average_jos_concentration("9",["2","10","20","40"])
    #check_average_jos_concentration("21",["15","20"])
    
    #check_average_sym_concentration("9",["2","10","20","40"])
    #check_average_sym_concentration("21",["15","20"])
    
    #check_dependance()
    #check_yue_dependance()    
    check_vanbaten_dependance()
    #give_overall_characteristics()     
     
    pylab.show()
