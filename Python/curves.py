#!/usr/bin/python
import numpy
import pylab
import math
import os

def Produce_Bunch():
    print os.getcwd()
    capillary_str=["3","5","8","10","20","40","60","80"]


    for dir_temp in capillary_str:
        dir_name="Results/"+dir_temp
        name=dir_name+"/concentration.dat"
        array=numpy.loadtxt(name)
        #print array
        new_conc=numpy.array(zip(array[1:,0]-array[:-1,0],array[1:,1]-array[:-1,1],array[1:,2]))
        print "Data=",new_conc
        coefficient=new_conc[:,1]/(new_conc[:,0]*(1.0-new_conc[:,2]))
        pylab.figure()
        pylab.plot(coefficient)


if __name__=="__main__":
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
    bubble_lengths=diameter*bubble_lengths_non
    slug_lengths=diameter*slug_lengths_non
    superficial_gas=bubble_velocities*holdups    
    print bubble_velocities/velocities_lbm    
    superficial_liquid=superficial_liq_lbm*bubble_velocities/velocities_lbm
    print "Liquid superficial=",superficial_liquid
    print "Bubble_velocities=",bubble_velocities
    print "Gas_superficial=",superficial_gas
    diffusion=1.118e-9    
    
    yue=2.0/hydraulic*numpy.sqrt(diffusion*bubble_velocities/(bubble_lengths+slug_lengths))*numpy.power(bubble_lengths/(bubble_lengths+slug_lengths),0.3)
    bercic=0.111*numpy.power(superficial_gas+superficial_liquid,1.19)/numpy.power((1-holdups)*(bubble_lengths+slug_lengths),0.57)    
    vanbaten=2.0/math.sqrt(math.pi)*numpy.sqrt(diffusion*bubble_velocities/(bubble_lengths-hydraulic))*4\
            *(bubble_lengths-hydraulic)/(hydraulic*(bubble_lengths+slug_lengths))\
            +2.0*math.sqrt(2)/math.pi*numpy.sqrt(diffusion*bubble_velocities/hydraulic)*4/(bubble_lengths+slug_lengths)
    
    fig=pylab.figure()    
    pylab.plot(bubble_velocities,yue,'k^',markersize=8)
    pylab.plot(bubble_velocities,bercic,'ko',markersize=8)
    pylab.plot(bubble_velocities,vanbaten,'ks',markersize=8)
    leg=["Yue","Bercic","Van Baten"]
    pylab.legend(leg)
    pylab.xlabel(r'''$U_{\mathrm{bubble}},\,\mathrm{m/s}$''',fontsize=30)
    pylab.ylabel(r'''$k_L a,\,\mathrm{s^{-1}}$''',fontsize=30)
    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    fig.subplots_adjust(left=0.17,bottom=0.17) 
    pylab.savefig("theoretical_correlations.eps",format="EPS",dpi=300)
    
    
    Produce_Bunch()
    pylab.show() 
