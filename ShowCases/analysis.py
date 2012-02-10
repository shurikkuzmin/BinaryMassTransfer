#!/usr/bin/python
import numpy
import os
import subprocess
import math
import pylab
import scipy
import scipy.optimize

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)


def give_overall_characteristics(flag,length_external):
    capillary_str=[str(x) for x in range(9,93,10)]
   
    #for length in range(200,1400,150):
    #    if not flag:
    if not flag:    
        length=length_external
        print length
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
            file_dir=str(length)+"/"+dir_name+"/"        
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
		    
            pylab.figure()
            pylab.imshow(phase)
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
    
    #capillaries=[]
    #reynolds=[]
    #velocities=[]
    #holdups=[]
    #bubble_lengths=[]
    #slug_lengths=[]
    #widths=[]
    #gas_velocities=[]
    #liq_velocities=[]
 
    #for counter,dir_name in enumerate(capillary_str):
    #    os.chdir("UnitCells/"+dir_name)
    #    geometry=numpy.loadtxt("geometry.dat").transpose()
    #    ux=numpy.loadtxt("vely0200000.dat").transpose()
    #    uy=numpy.loadtxt("velx0200000.dat").transpose()
    #    dims=geometry.shape
       
    #    center=geometry[dims[0]/2,:]
    #    z1 = numpy.min(numpy.where(center < 0.0))
    #    z2 = numpy.max(numpy.where(center < 0.0))
    #    if z1==0:
    #        z2=numpy.min(numpy.where(center>0.0))+dims[1]
    #        z1=numpy.max(numpy.where(center>0.0))
    #   
    #    holdups.append(len(numpy.where(geometry<0)[0])/float((dims[0]-2)*dims[1]))
    #    velocities.append(ux[0,0])
    #    liq_velocities.append(ux[0,0]-numpy.sum(ux[1:-1,((z2+z1+dims[1])/2)%dims[1]])/float(dims[0]-2))
    #    gas_velocities.append(velocities[-1]*holdups[-1])

    #    prof=geometry[:,((z1+z2)/2)%dims[1]]
    #    widths.append(Get_Zero(prof))     
    #    bubble_lengths.append((z2-z1)/float(dims[0]-2))
    #    slug_lengths.append(15.0-(z2-z1)/float(dims[0]-2))
    #    capillaries.append(velocities[-1]*(2.0/3.0)/math.sqrt(8*0.04*0.04/9.0))
    #    reynolds.append(velocities[-1]*(dims[0]-2)/(2.0/3.0))

     
    #    os.chdir("../..")
    
    #print "Bubble velocities=",velocities
    #print "Liquid velocities=",liq_velocities
    #print "Gas velocities=",gas_velocities
    #print "Holdups=",holdups
    #print "Widths=",widths
    #print "Capillaries=",capillaries
    #print "Reynolds=",reynolds_orig
    #print "Bubble lengths=",bubble_lengths
    #print "Slug lengths=",slug_lengths

def check_forcing():
    capillaries200= [0.073722338265551551, 0.1719263011747692, 0.28077964694645074, 0.38356059451146424,\
                     0.49511917250675719, 0.60226447292781515, 0.71706884817900307, 0.83144814548157386,\
                     0.9389104954359766]
    capillaries350= [0.081920737190919246, 0.19092818276110241, 0.30719805520107568, 0.4301302757125734,\
                     0.55270779219396804, 0.67283957781686643, 0.79980375053962127, 0.92506757717322885,\
                     1.0487547283956218]
    capillaries500= [0.093085865655324646, 0.21529256443354303, 0.34686116743702095, 0.48753362635204572,\
                     0.6240645822396409, 0.75835677866103002, 0.89734106273915526, 1.1484222582865664,\
                     1.4197024992681524]
    capillaries650= [0.10376071224615181, 0.24899292714257107, 0.39725308818088317, 0.55338978802236449,\
                     0.71158813458225767, 0.87306468227340051, 1.933969836854418, 2.702369675048951, 1.7669263331111857]
    capillaries800= [0.12101798570494038, 0.28984689346506176, 0.47031523593501479, 0.65124792517512475,\
                     0.83306997849960196, 1.7364136585066581, 1.5713088393249, 2.7530947436726914, 2.8255210307803655]
    capillaries950= [0.14337631171472745, 0.34805495902191214, 0.55946628451387026, 0.77332791142578172,\
                     0.99267743490189853, 1.6104939360630546, 2.207150958701189, 2.4290665048039806, 3.6098015344539602]
    capillaries1100= [0.1729303453571891, 0.43306952220305733, 0.69650644938457806, 0.95996530041184946,\
                      1.6575726302019211, 2.1152879307565264, 2.6225967942043553, 3.1339242816077229, 3.9810722367508178]

    capillaries=numpy.array(zip(capillaries200,capillaries350,capillaries500,capillaries650,capillaries800,capillaries950,capillaries1100))
    length=200+150*numpy.arange(0,7)
    print capillaries
    pylab.plot(length,capillaries[0,:])
    pylab.plot(length,capillaries[1,:])
    pylab.plot(length,capillaries[2,:])
    pylab.plot(length,capillaries[3,:])
    pylab.plot(length,capillaries[4,:])
    
    x=numpy.arange(0,7)
    pylab.plot(length,0.4951*(1.0+(length-200.0)*(length-200.0)/(150*150*26.0)),"+")
if __name__=="__main__":

    #give_overall_characteristics(False,1100)     
    check_forcing() 
     
    pylab.show()
