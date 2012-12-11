import os
import numpy
import pylab
import glob
def get_scales(dir_name):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    scales=[0.3,0.5,1,2,4,6,8,10,15,20]    
    scale_str=["03","05","1","2","4","6","8","10","15","20"]
    flag=False
    for scale_counter,scale in enumerate(scales):
        if velocities[counter]*scale>0.12:
            scale3=scales[scale_counter-1]
            scale2=scales[scale_counter-2]
            scale1=scales[scale_counter-3]
            flag=True 
            break
    if not flag:
        scale1=scales[-3]
        scale2=scales[-2]
        scale3=scales[-1] 
    return scale1,scale2,scale3 
    
    
def copy_one_directory(dir_name,scale):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    units=[4,6,8,10]
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    scales=[0.3,0.5,1,2,4,6,8,10,15,20]    
    scale_str=["03","05","1","2","4","6","8","10","15","20"]
    counter=capillary_str.index(dir_name)

    os.chdir(dir_name)    
    os.chdir(str(scale))
    for unit in units:
        os.chdir(str(unit))            
        subprocess.call(['scp','shurik@checkers.westgrid.ca:/home/shurik/UnitCells/'+dir_name+"/"+dir_name_scale+"/"+str(unit)+"/density*.dat","."])
        print "Done with ",unit
        os.chdir("..")
    os.chdir("..")
    print "Scale ",scale," is finished"
def show_one_image(dir_name,scale,unit,time_step):
    os.chdir("UnitCells/"+dir_name)    
    os.chdir(str(scale))
    os.chdir(str(unit))            
    arr=numpy.loadtxt("density"+time_step+".dat")        
    pylab.imshow(arr)
    pylab.colorbar()

def analyze_one_unit(cluster_name,dir_name,scale,unit,main_flag):
    shift=0
    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    capillary_str=["9","21","42","60","84"]
    holdups=[0.30,0.28,0.27,0.25,0.22]
    velocities_liq=numpy.array([ 0.00459266431535,0.0107317293541,0.0208964048893,0.0292353942998,0.0362306511905])    
    velocities_gas=holdups*velocities    

    
    if cluster_name=="bugaboo":
        os.chdir("UnitCells/Bugaboo/"+dir_name)
        ux=numpy.loadtxt("../../"+dir_name+"/vely0200000.dat").transpose()
        uy=numpy.loadtxt("../../"+dir_name+"/velx0200000.dat").transpose()
        geometry=numpy.loadtxt("../../"+dir_name+"/geometry.dat").transpose()

    else:
        os.chdir("UnitCells/"+dir_name)
        ux=numpy.loadtxt("vely0200000.dat").transpose()
        uy=numpy.loadtxt("velx0200000.dat").transpose()
        geometry=numpy.loadtxt("geometry.dat").transpose()

        
    ind_liquid=numpy.where(geometry>0)
    os.chdir(str(scale))
    os.chdir(str(unit))            
    files_names=glob.glob("density*.dat")
    files_names=sorted(files_names)
    concentrations=numpy.zeros([len(files_names),unit+1])
    concentrations_real=numpy.zeros([len(files_names),unit+1])
    concentrations_center=numpy.zeros([len(files_names),unit+1])
    average_concentrations=numpy.zeros([len(files_names),unit])
    print concentrations.shape    
    if main_flag:    
        for time_counter,file_name in enumerate(files_names):
            arr=numpy.loadtxt(file_name)
            dims=arr.shape
            for un in range(0,unit):
                concentrations[time_counter,un]=numpy.sum(arr[:,un*3000+shift]*ux[:,shift])/numpy.sum(ux[:,shift])
                concentrations_real[time_counter,un]=numpy.sum(arr[:,un*3000+shift])/200
                concentrations_center[time_counter,un]=arr[100,un*3000+shift]
                cross=numpy.array(ind_liquid)
                print cross.shape
                cross[1,:]=3000*un+cross[1,:]                
                average_concentrations[time_counter,un]=numpy.sum(arr[cross.tolist()])/len(ind_liquid[0])
            
            concentrations[time_counter,-1]=numpy.sum(arr[:,-1]*ux[:,0])/numpy.sum(ux[:,0])
            concentrations_real[time_counter,-1]=numpy.sum(arr[:,-1])/200
            concentrations_center[time_counter,-1]=arr[100,-1]
            print file_name
        numpy.savetxt("concentrations.dat",concentrations)
        numpy.savetxt("concentrations_real.dat",concentrations_real)
        numpy.savetxt("concentrations_center.dat",concentrations_center)
        numpy.savetxt("average_concentrations.dat",average_concentrations)
    else:
        concentrations=numpy.loadtxt("concentrations.dat")
        concentrations_real=numpy.loadtxt("concentrations_real.dat")
        concentrations_center=numpy.loadtxt("concentrations_center.dat")
        average_concentrations=numpy.loadtxt("average_concentrations.dat")        
        
    legs=[str(x) for x in range(0,unit+1)]
    legs_segments=[str(x+2)+" segment" for x in range(0,unit-2)]
    for un in range(0,unit+1):
        pylab.figure(1)
        pylab.plot(concentrations_real[:,un])
        pylab.figure(2)        
        pylab.plot(concentrations[:,un])
        pylab.figure(3)
        pylab.plot(concentrations_center[:,un])
    pylab.figure(1)
    pylab.legend(legs)
    pylab.figure(2)
    pylab.legend(legs)
    pylab.figure(3)
    pylab.legend(legs)
    

    pylab.figure(4)
    coefficients=numpy.zeros([len(files_names),unit-2])
    all_times=50000*numpy.arange(0,len(files_names))
    counter=capillary_str.index(dir_name)
    for un in range(0,unit-2):
        coefficients[:,un]=3000.0/(all_times*velocities[counter])*numpy.log((1-numpy.abs(concentrations[:,un+1]))/(1-numpy.abs(concentrations[:,un+2])))
        pylab.plot(coefficients[:,un])
    #pylab.colorbar()
    pylab.legend(legs_segments)
    
    pylab.figure(5)
    for un in range(0,unit-2):
        pylab.plot(numpy.log((1-numpy.abs(concentrations_center[:,un+1]))/(1-numpy.abs(concentrations_center[:,un+2]))))
    #pylab.colorbar()
    pylab.legend(legs_segments)
    
    fig=pylab.figure(6)
    styles=["k+","k.","ko","ks","kv","k^","kd","k--","kD","kd","kx"]
    for un in range(0,unit):
        pylab.plot(all_times/10000,average_concentrations[:,un],styles[un])
    legs=[str(x)+" unit" for x in range(0,unit)]    
    pylab.legend(legs,loc=2)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$C_{aver}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale}='''+str(scale)+'''$''',fontsize=30)
    fig.subplots_adjust(left=0.15,bottom=0.15)
    fig.savefig("aver_units"+str(unit)+"scale"+str(scale)+".eps",format="EPS",dpi=300)
    
    fig=pylab.figure(7)
    styles=["k+","k.","ko","ks","kv","k^","kd","k--","kD","kd","kx"]
    for un in range(0,unit-1):
        pylab.plot(all_times/10000,numpy.log((1-numpy.abs(concentrations[:,un]))/(1-numpy.abs(concentrations[:,un+1])))
,styles[un])
    legs=[str(x)+"-"+str(x+1)+" segment" for x in range(0,unit-1)]    
    pylab.legend(legs,loc=2)
    
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$k_L a\frac{L}{U}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale}='''+str(scale)+'''$''',fontsize=30)
    fig.subplots_adjust(left=0.2,bottom=0.15)
    fig.savefig("coeff_units"+str(unit)+"scale"+str(scale)+".eps",format="EPS",dpi=300)
  
    
    pylab.figure(8)
    pylab.plot(3000.0/(velocities[counter]*scale*all_times)*numpy.log(1.0/(1.0-average_concentrations[:,0])))
    pylab.plot(3000.0/(velocities[counter]*scale*all_times)*numpy.log(1.0/(1.0-average_concentrations[:,1])))
    #pylab.ylim(ymax=0.3,ymin=0.0)
    
    #This section is to use a new definition
    fig=pylab.figure(9)
    window_len=2
    
    mass=average_concentrations*3000.0*200.0/(1.0-holdups[counter])
    inlet_fluxes=numpy.zeros_like(average_concentrations) 
    outlet_fluxes=numpy.zeros_like(average_concentrations)
    for un in range(0,unit):    
        inlet_fluxes[:,un]=concentrations[:,un]*numpy.sum(ux[1:-1,0])*scale
        outlet_fluxes[:,un]=concentrations[:,un+1]*numpy.sum(ux[1:-1,0])*scale
    print inlet_fluxes.shape,outlet_fluxes.shape
    print inlet_fluxes[window_len:,:].shape,inlet_fluxes[:-window_len,:].shape
    print outlet_fluxes[window_len:,:].shape,outlet_fluxes[:-window_len,:].shape

    outcome=numpy.zeros([len(all_times)-window_len,unit])
    for un in range(0,unit):
        outcome[:,un]=((mass[window_len:,un]-mass[:-window_len,un])/(50000*window_len)\
                 +inlet_fluxes[window_len/2:-window_len/2,un]\
                 -outlet_fluxes[window_len/2:-window_len/2,un])\
            /(200.0*scale*(velocities_liq[counter]+velocities_gas[counter])*(1.0-average_concentrations[window_len/2:-window_len/2,un]))
        pylab.plot(all_times[window_len/2:-window_len/2]/10000,outcome[:,un],styles[un])
    pylab.ylim(ymin=0.0,ymax=0.5)
    legs=[str(x)+" unit" for x in range(0,unit)]    
   
    pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$k_L a\frac{L}{U_{\mathrm{gas}}+U_{\mathrm{liq}}}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale}='''+str(scale)+'''$''',fontsize=30)
    fig.subplots_adjust(left=0.2,bottom=0.15)
    fig.savefig("right_def_"+str(unit)+"scale"+str(scale)+".eps",format="EPS",dpi=300)
    
    
    #old_concentrations=numpy.loadtxt("concentration.dat")
    #flux[:,un]=(concentrations[:,un+1]-concentrations[:,un])*numpy.sum(ux[1:-1,0])\
    #        /(200.0*(velocities_liq[counter]+velocities_gas[counter])*(1.0-average_concentrations[:,un]))
    


def analyze_increasing_diffusion_unit(dir_name,scale,scale_diffusion,unit,main_flag,start_iteration,last_iteration):
    shift=0
    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    capillary_str=["9","21","42","60","84"]
    os.chdir("IncreasingDiffusion/"+dir_name)
    ux=numpy.loadtxt("vely0200000.dat").transpose()
    uy=numpy.loadtxt("velx0200000.dat").transpose()
    geometry=numpy.loadtxt("geometry.dat").transpose()
       
    ind_liquid=numpy.where(geometry>0)
    os.chdir(str(scale))
    os.chdir(str(scale_diffusion))
    os.chdir(str(unit))            
    files_names=glob.glob("density*.dat")
    files_names=sorted(files_names)
    #concentrations=numpy.zeros([len(files_names),unit+1])
    #concentrations_real=numpy.zeros([len(files_names),unit+1])
    #concentrations_center=numpy.zeros([len(files_names),unit+1])
    #average_concentrations=numpy.zeros([len(files_names),unit])
    concentrations=numpy.zeros([last_iteration,unit+1])
    concentrations_real=numpy.zeros([last_iteration,unit+1])
    concentrations_center=numpy.zeros([last_iteration,unit+1])
    average_concentrations=numpy.zeros([last_iteration,unit])
    
    print concentrations.shape    
    if main_flag:    
        for time_counter,file_name in enumerate(files_names[start_iteration:last_iteration]):
            time_counter=time_counter+start_iteration            
            arr=numpy.loadtxt(file_name)
            dims=arr.shape
            for un in range(0,unit):
                concentrations[time_counter,un]=numpy.sum(arr[:,un*3000+shift]*ux[:,shift])/numpy.sum(ux[:,shift])
                concentrations_real[time_counter,un]=numpy.sum(arr[:,un*3000+shift])/200
                concentrations_center[time_counter,un]=arr[100,un*3000+shift]
                cross=numpy.array(ind_liquid)
                print cross.shape
                cross[1,:]=3000*un+cross[1,:]                
                average_concentrations[time_counter,un]=numpy.sum(arr[cross.tolist()])/len(ind_liquid[0])
            
            concentrations[time_counter,-1]=numpy.sum(arr[:,-1]*ux[:,0])/numpy.sum(ux[:,0])
            concentrations_real[time_counter,-1]=numpy.sum(arr[:,-1])/200
            concentrations_center[time_counter,-1]=arr[100,-1]
            print file_name
        numpy.savetxt("concentrations.dat",concentrations)
        numpy.savetxt("concentrations_real.dat",concentrations_real)
        numpy.savetxt("concentrations_center.dat",concentrations_center)
        numpy.savetxt("average_concentrations.dat",average_concentrations)
    else:
        concentrations=numpy.loadtxt("concentrations.dat")
        concentrations_real=numpy.loadtxt("concentrations_real.dat")
        concentrations_center=numpy.loadtxt("concentrations_center.dat")
        average_concentrations=numpy.loadtxt("average_concentrations.dat")        
        
    legs=[str(x) for x in range(0,unit+1)]
    legs_segments=[str(x+2)+" segment" for x in range(0,unit-2)]
    #for un in range(0,unit+1):
        #pylab.figure(1)
        #pylab.plot(concentrations_real[:,un])
        #pylab.figure(2)        
        #pylab.plot(concentrations[:,un])
        #pylab.figure(3)
        #pylab.plot(concentrations_center[:,un])
    #pylab.figure(1)
    #pylab.legend(legs)
    #pylab.figure(2)
    #pylab.legend(legs)
    #pylab.figure(3)
    #pylab.legend(legs)
    

    pylab.figure(4)
    coefficients=numpy.zeros([last_iteration,unit-2])
    #all_times=50000*numpy.arange(0,len(files_names))
    all_times=50000*numpy.arange(0,last_iteration)
    print concentrations.shape
    print all_times.shape
    counter=capillary_str.index(dir_name)
    for un in range(0,unit-2):
        coefficients[:,un]=3000.0/(all_times*velocities[counter])*numpy.log((1-numpy.abs(concentrations[:,un+1]))/(1-numpy.abs(concentrations[:,un+2])))
        pylab.plot(coefficients[:,un])
    #pylab.colorbar()
    pylab.legend(legs_segments)
    
    pylab.figure(5)
    for un in range(0,unit-2):
        pylab.plot(numpy.log((1-numpy.abs(concentrations_center[:,un+1]))/(1-numpy.abs(concentrations_center[:,un+2]))))
    #pylab.colorbar()
    pylab.legend(legs_segments)
    
    fig=pylab.figure(6)
    styles=["k+","k.","ko","ks","kv","k^","kd","k--","kD","kd","kx"]
    for un in range(0,unit):
        pylab.plot(all_times/10000,average_concentrations[:,un],styles[un])
    legs=[str(x)+" unit" for x in range(0,unit)]    
    pylab.legend(legs,loc=2)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$C_{aver}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale_U}='''+str(scale)+''',\,\mathrm{Scale_D}='''+str(scale*scale_diffusion)+r'''$''',fontsize=30)
    fig.subplots_adjust(left=0.15,bottom=0.15)
    fig.savefig("aver_units"+str(unit)+"scaleu"+str(scale)+"scaled"+str(scale_diffusion)+".eps",format="EPS",dpi=300)
    
    fig=pylab.figure(7)
    styles=["k+","k.","ko","ks","kv","k^","kd","k--","kD","kd","kx"]
    for un in range(0,unit-1):
        pylab.plot(all_times/10000,numpy.log((1-numpy.abs(concentrations[:,un]))/(1-numpy.abs(concentrations[:,un+1])))
,styles[un])
    legs=[str(x)+"-"+str(x+1)+" segment" for x in range(0,unit-1)]    
    pylab.legend(legs,loc=2)
    
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$k_L a\frac{L}{U}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale_U}='''+str(scale)+''',\,\mathrm{Scale_D}='''+str(scale*scale_diffusion)+r'''$''',fontsize=30)
    fig.subplots_adjust(left=0.2,bottom=0.15)
    fig.savefig("novortex"+str(unit)+"scaleu"+str(scale)+"scaled"+str(scale_diffusion)+".eps",format="EPS",dpi=300)
  
    print "Peclet=",scale*velocities[counter]*200.0/(1.0/3.0*(1.0/1.99-0.5)*scale*scale_diffusion)
    #pylab.figure(8)
    #pylab.plot(3000.0/(velocities[counter]*scale*all_times)*numpy.log(1.0/(1.0-average_concentrations[:,0])))
    #pylab.plot(3000.0/(velocities[counter]*scale*all_times)*numpy.log(1.0/(1.0-average_concentrations[:,1])))

def analyze_flux(dir_name,scale,scale_diffusion,unit,main_flag,start_iteration,last_iteration):
    holdups=numpy.array([0.302466666667,0.288423333333,0.271146666667,0.252866666667,0.229446666667])
    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    velocities_liq=numpy.array([ 0.00459266431535,0.0107317293541,0.0208964048893,0.0292353942998,0.0362306511905])    
    velocities_gas=holdups*velocities    
    capillary_str=["9","21","42","60","84"]
    styles=["k+","k.","ko","ks","kv","k^","kd","k--","kD","kd","kx"]
    os.chdir("IncreasingDiffusion/"+dir_name)
    ux=numpy.loadtxt("vely0200000.dat").transpose()
    uy=numpy.loadtxt("velx0200000.dat").transpose()
    geometry=numpy.loadtxt("geometry.dat").transpose()
    counter=capillary_str.index(dir_name)   
    ind_liquid=numpy.where(geometry>0)
    os.chdir(str(scale))
    os.chdir(str(scale_diffusion))
    os.chdir(str(unit))            
    #concentrations=numpy.zeros([last_iteration,unit+1])
    #concentrations_real=numpy.zeros([last_iteration,unit+1])
    #concentrations_center=numpy.zeros([last_iteration,unit+1])
    #average_concentrations=numpy.zeros([last_iteration,unit])
    concentrations=numpy.loadtxt("concentrations.dat")
    concentrations_real=numpy.loadtxt("concentrations_real.dat")
    concentrations_center=numpy.loadtxt("concentrations_center.dat")
    average_concentrations=numpy.loadtxt("average_concentrations.dat")        
    
    flux=numpy.zeros_like(concentrations)
    flux_aver=numpy.zeros_like(concentrations[:,:-1])
    print concentrations.shape,average_concentrations.shape
    for un in range(0,unit):
        flux[:,un]=(concentrations[:,un+1]-concentrations[:,un])*numpy.sum(ux[1:-1,0])\
            /(200.0*(velocities_liq[counter]+velocities_gas[counter])*(1.0-average_concentrations[:,un]))
    pylab.figure(1)
    for un in range(0,unit+1):
        pylab.plot(concentrations[:,un])
    pylab.figure(2)
    for un in range(0,unit+1):
        pylab.plot(flux[:,un],styles[un])
    pylab.figure(3)
    ux_aver=numpy.sum(ux[numpy.where(geometry>0)])/len(numpy.where(geometry>0)[0])
    for un in range(0,unit-1):
        flux_aver[:,un]=(average_concentrations[:,un+1]-average_concentrations[:,un])*ux_aver\
            /((velocities_liq[counter]+velocities_gas[counter])*(1.0-average_concentrations[:,un]))
    pylab.plot(flux_aver)
    
    fig=pylab.figure(2)
    legs=[str(x)+" unit" for x in range(0,unit)]    
    pylab.legend(legs,loc=1)
        
    pylab.ylim(ymax=0.5,ymin=0.0)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
 
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$k_L a \frac{L}{U_{\mathrm{gas}}+U_{\mathrm{liq}}}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale_U}='''+str(scale)+''',\, Scale_D='''+str(scale_diffusion)+'''$''',fontsize=30)
    fig.subplots_adjust(left=0.2,bottom=0.15)
    fig.savefig("flux_moving_window"+str(unit)+"scaleu"+str(scale)+"scaled"+str(scale_diffusion)+".eps",format="EPS",dpi=300)

    
def moving_window(cluster_name,dir_name,scale,unit):
    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    capillary_str=["9","21","42","60","84"]
    counter=capillary_str.index(dir_name)
    window_len=2

    if cluster_name=="bugaboo":
        os.chdir("UnitCells/Bugaboo/"+dir_name+"/"+str(scale)+"/"+str(unit))
    else:
        os.chdir("UnitCells/"+dir_name+"/"+str(scale)+"/"+str(unit))

    average_concentrations=numpy.loadtxt("average_concentrations.dat")        
    dims=average_concentrations.shape    
    files_names=glob.glob("density*.dat")
    all_times=50000*numpy.arange(0,average_concentrations.shape[0])
    x=average_concentrations[:,0]    
    #pylab.plot(average_concentrations)
    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    w=numpy.ones(window_len,'d')
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    
    styles=["k+","k.","ko","ks","kv","k^","kd","k--","kD","kd","kx"]
    
    fig=pylab.figure()
    for un in range(0,dims[1]):
        conc=[]
        for i in range(0,dims[0]-window_len):
            conc.append(3000.0/(velocities[counter]*scale*(50000*window_len))\
                        *numpy.log((1.0-average_concentrations[i,un])/(1.0-average_concentrations[i+window_len,un])))   
        pylab.plot(conc,styles[un])
 
    legs=[str(x)+" unit" for x in range(0,unit)]    
    pylab.legend(legs,loc=1)
        
    pylab.ylim(ymax=0.5,ymin=0.0)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
 
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$k_L a \frac{L}{U_{\mathrm{gas}}+U_{\mathrm{liq}}}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale}='''+str(scale)+'''$''',fontsize=30)
    fig.subplots_adjust(left=0.2,bottom=0.15)
    fig.savefig("aver_moving_window"+str(unit)+"scale"+str(scale)+".eps",format="EPS",dpi=300)
 
def moving_increasing_diffusion_window(dir_name,scale,scale_diffusion,unit):
    velocities=numpy.array([0.0055,0.0143,0.0297,0.0424,0.05538])
    capillary_str=["9","21","42","60","84"]
    counter=capillary_str.index(dir_name)
    window_len=2

    os.chdir("IncreasingDiffusion/"+dir_name+"/"+str(scale)+"/"+str(scale_diffusion)+"/"+str(unit))
    
    average_concentrations=numpy.loadtxt("average_concentrations.dat")        
    dims=average_concentrations.shape    
    files_names=glob.glob("density*.dat")
    #all_times=50000*numpy.arange(0,average_concentrations.shape[0])
    #x=average_concentrations[:,0]    
    #pylab.plot(average_concentrations)
    #s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #w=numpy.ones(window_len,'d')
    #y=numpy.convolve(w/w.sum(),s,mode='valid')
    
    styles=["k+","k.","ko","ks","kv","k^","kd","k--","kD","kd","kx"]
    
    fig=pylab.figure()
    for un in range(0,dims[1]):
        conc=[]
        for i in range(0,dims[0]-window_len):
            conc.append(3000.0/(velocities[counter]*scale*(50000*window_len))\
                        *numpy.log((1.0-average_concentrations[i,un])/(1.0-average_concentrations[i+window_len,un])))   
        pylab.plot(conc,styles[un])
 
    legs=[str(x)+" unit" for x in range(0,unit)]    
    pylab.legend(legs,loc=1)
        
    pylab.ylim(ymax=0.5,ymin=0.0)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
 
    pylab.xlabel(r'''$\mathrm{Time}\times 10^4$''',fontsize=30)
    pylab.ylabel(r'''$k_L a \frac{L}{U_{\mathrm{gas}}+U_{\mathrm{liq}}}$''',fontsize=30)    
    pylab.title(r'''$\mathrm{Scale_U}='''+str(scale)+''',\, Scale_D='''+str(scale_diffusion)+'''$''',fontsize=30)
    fig.subplots_adjust(left=0.2,bottom=0.15)
    fig.savefig("aver_moving_window"+str(unit)+"scaleu"+str(scale)+"scaled"+str(scale_diffusion)+".eps",format="EPS",dpi=300)


def analyze_average_diffusion(dir_name,scale,scale_diffusion):
    capillary_str=[str(x) for x in [9,21,42,60,84]]
    capillaries=[0.097,0.254,0.526,0.750,1.040]
    velocities=[0.0055,0.0143,0.0297,0.0424,0.05538]
    holdups=[0.30,0.28,0.27,0.25,0.22]
    #scale=[8,10,20,40]
    os.chdir("IncreasingDiffusion/"+dir_name+"/"+str(scale)+"/"+str(scale_diffusion))
    
    counter=capillary_str.index(dir_name) 
    capillary=capillaries[counter]
    gas_holdup=holdups[counter]
    #styles=["k+","ko","ks","k^","kv","kd"]    
    
    styles=['k-','k--','k-.','k:','k.','k,']    
    
    legs=[]
    concentration=numpy.loadtxt("concentration.dat")
    aver_concentration=concentration[1:,1]/(3000.0*200.0*(1.0-gas_holdup))
    pylab.figure(99)
    pylab.plot(aver_concentration) #,styles[scale_counter])
    
    fig=pylab.figure(100)
        
    ##this plot is in terms of unit cells        
    ##pylab.plot(velocities[counter]*scale[scale_counter]*concentration[1:,0]/3000.0,3000.0/(velocities[counter]*scale[scale_counter]*concentration[1:,0])*numpy.log(1.0/(1.0-aver_concentration)),styles[scale_counter],linewidth=2)                
    ##this plot is in terms of time        
    pylab.plot(concentration[1:,0]/1000.0,3000.0/(velocities[counter]*scale*concentration[1:,0])*numpy.log(1.0/(1.0-aver_concentration)),"ks")                
    fig.subplots_adjust(top=0.9,bottom=0.15,left=0.18)
    ##pylab.legend(legs)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.ylim(ymax=0.5,ymin=0.0)
    ##pylab.xlabel(r'''$\mathrm{Cell\,Units}$''',fontsize=30)
    pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    
    pylab.ylabel(r'''$k_L a \frac{L_{\mathrm{unit}}}{U_{\mathrm{liq}}+U_{\mathrm{gas}}}$''',fontsize=30)    
    pylab.title(r'''$Ca='''+str(capillary)+''',\,'''+"Pe=2644"+'''$''',fontsize=30)
    pylab.savefig("volume_ca_"+str(capillary)+"_scaleu"+str(scale)+"scaled"+str(scale_diffusion)+".eps",format="EPS",dpi=300)
    #pylab.figure(101)
    #window_len=10     
    #conc=[]
    #for i in range(0,len(aver_concentration)-window_len):
        #conc.append(3000.0/(velocities[counter]*scale[scale_counter]*1000.0*window_len)\
                    #*numpy.log((1.0-aver_concentration[i])/(1.0-aver_concentration[i+window_len])))   
    #pylab.plot(conc) #styles[scale_counter])
    #pylab.ylim(ymin=0.0,ymax=0.5)
    #pylab.figure(102)
    #vanbaten=[]
    #mult_vanbaten=velocities[counter]*scale
    #for i in range(0,len(aver_concentration)-window_len):
        #vanbaten.append((aver_concentration[i+window_len]-aver_concentration[i])*3000.0/(1000*mult_vanbaten*window_len*(1.0-aver_concentration[i+window_len/2])))   
    #pylab.plot(vanbaten) #,styles[scale_counter])
    #pylab.ylim(ymin=0.0,ymax=0.5)
    ##pylab.figure(99,figsize=(12,3))
    #print 3000.0/(velocities[counter]*scale[scale_counter]*concentration[exam_point,0])*numpy.log(1.0/(1.0-aver_concentration[exam_point-1]))
    #legs.append(scale_str[scale_counter])        
    ##cont=pylab.contour(arr,levels=c_levels,linestyles="dotted",colors=['k'],linewidths=[2],extent=(0,15,0,1))
    ##if scale_counter==len(scale_str)-subtract-1:        
    ##    pylab.clabel(cont)
    ##else:
    ##    cont.clabel(fmt=scale_str[scale_counter],fontsize="12",inline_spacing=10)
        
    ##fig=pylab.figure(99)
    
    #fig=pylab.figure(100)    
    #fig.subplots_adjust(top=0.9,bottom=0.15,left=0.18)
    ##pylab.legend(legs)
    #pylab.xticks(fontsize=20)
    #pylab.yticks(fontsize=20)
    #pylab.ylim(ymax=0.5,ymin=0.0)
    ##pylab.xlabel(r'''$\mathrm{Cell\,Units}$''',fontsize=30)
    #pylab.xlabel(r'''$\mathrm{Time}\times 10^3$''',fontsize=30)
    
    #pylab.ylabel(r'''$k_L a \frac{L_{\mathrm{unit}}}{U_{\mathrm{liq}}+U_{\mathrm{gas}}}$''',fontsize=30)    
    #pylab.title(r'''$Ca='''+str(capillary)+'''$''',fontsize=30)
    ##pylab.savefig("aver_conc_scale_ca"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)
    ##pylab.savefig("aver_conc_scale_ca_time"+str(capillary)[0:1]+str(capillary)[3:]+".eps",format="EPS",dpi=300)


    
if __name__=="__main__":
    import os.path
    curdir=os.getcwd()
    print curdir
    
    #analyze_one_scale("9",10)
    #show_one_image("9",10,4,"0600000")      
    #analyze_one_unit("bugaboo","9",20,4,False)
    analyze_one_unit("bugaboo","9",10,6,False)    
    
    #analyze_increasing_diffusion_unit("84",2,5,4,False,0,600000/50000)
        
    #analyze_increasing_diffusion_unit("84",2,5,6,False,0,2000000/50000)
    #analyze_flux("84",2,5,6,False,0,2000000/50000)    
    #os.chdir(curdir)
    #analyze_average_diffusion("84",2,5)
    #os.chdir(curdir)    
    #moving_increasing_diffusion_window("84",2,5,6)

    #analyze_increasing_diffusion_unit("84",2,5,8,True,0,1350000/50000)
    
    #analyze_increasing_diffusion_unit("84",2,10,4,False,0,1250000/50000)
    #os.chdir(curdir)    
    #analyze_average_diffusion("84",2,10)
    #analyze_flux("84",2,10,6,True,100000/50000,2000000/50000)
    #os.chdir(curdir)    
    
    #analyze_increasing_diffusion_unit("84",2,5,6,False,100000/50000,2000000/50000)
    #analyze_average_diffusion("84",2,10)    
    #analyze_increasing_diffusion_unit("84",2,10,6,True,100000/50000,2000000/50000)
    #analyze_increasing_diffusion_unit("84",2,10,8,True,50000/50000,2000000/50000)
    
    #analyze_increasing_diffusion_unit("84",2,20,4,True,100000/50000,1250000/50000)
    #analyze_increasing_diffusion_unit("84",2,20,6,True,50000/50000,2000000/50000) - corrupted
    #analyze_increasing_diffusion_unit("84",2,20,8,True,50000/50000,2000000/50000) - corrupted

    #analyze_increasing_diffusion_unit("84",2,40,6,False,0,2000000/50000)
    #os.chdir(curdir)
    #analyze_average_diffusion("84",2,40)
    #os.chdir(curdir)
    #analyze_increasing_diffusion_unit("84",2,40,6,False,0,2000000/50000)
    #analyze_increasing_diffusion_unit("84",2,20,8,True,50000/50000,2000000/50000) - corrupted
    

    
    #moving_increasing_diffusion_window("84",2,5,6)
    #moving_increasing_diffusion_window("84",2,5,8)
    
    
    #moving_window("bugaboo","9",20,6)
    pylab.show()              
