//Jos mass transfer
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>

//Domain size
int NY;
int NX;
int NUM;
int radius;

//Other constants
const int NPOP=9;

//Time steps
int N=5000;
int NOUTPUT=100;

//Fields and populations
double *f;
double *f2;
double *rho;
double *ux;
double *uy;
int * geometry;

//Boundary conditions
double conc_wall=1.0;

//BGK relaxation parameter
double omega=1.0;
double omega_high=1.5; //0.2; ///1.99;
double omega_low=1.0; //1.99;
double omega_plus=2.0-omega;
double omega_minus=omega;


//Diffusion parameter
double ce=1.0/3.0;


//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
double weights_trt[]={0.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};
int cx[]={0,1,0,-1,0,1,-1,-1,1};
int cy[]={0,0,1,0,-1,1,1,-1,-1};
int compliment[]={0,3,4,1,2,7,8,5,6};
int pxx[]={0, 1, -1, 1, -1, 0, 0, 0, 0};
int pxy[]={0, 0, 0, 0, 0, 1, -1, 1, -1};
void writedensity(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iX=0; iX<NX; ++iX)
	{
		for (int iY=0; iY<NY; iY++)
		{
			int counter=iY*NX+iX;
			fout<<rho[counter]<<" ";
		}
		fout<<"\n";
	}
}

void writegeometry(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iX=0; iX<NX; ++iX)
	{
		for (int iY=0; iY<NY; iY++)
		{
			int counter=iY*NX+iX;
			fout<<geometry[counter]<<" ";
		}
		fout<<"\n";
	}
}



void init()
{
	//Creating arrays
	f=new double[NUM*NPOP];
	f2=new double[NUM*NPOP];
	
	//Bulk nodes initialization
	double feq;
	double geq;
	double sum;
	
	for(int iY=0;iY<NY;iY++)
		for(int iX=0;iX<NX;iX++)
		{
			int  counter=iY*NX+iX;
			double dense_temp=rho[counter];
			double ux_temp=ux[counter];
			double uy_temp=uy[counter];
            
            for (int k=0; k<NPOP; k++)
			{
				feq=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
								                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
								                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				f[counter*NPOP+k]=feq;
			}
			
		}

}

void collide_column_bgk(int coor_y,int coor_bottom,int coor_top)
{
	for(int iX=coor_bottom;iX<=coor_top;iX++)
	{
		int counter=coor_y*NX+iX;
		rho[counter]=0.0;

		double sum=0;
		for(int k=0;k<NPOP;k++)
			sum+=f[counter*NPOP+k];
	
		rho[counter]=sum;

		double dense_temp=rho[counter];
		double ux_temp=ux[counter];
		double uy_temp=uy[counter];

		double feqeq[NPOP];

		for (int k=0; k<NPOP; k++)
			feqeq[k]=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
            						         +4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
															+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
					    			                         +2.0*ux_temp*uy_temp*cx[k]*cy[k]));

		for(int k=0; k < NPOP; k++)
			f2[counter*NPOP+k]=f[counter*NPOP+k]*(1.0-omega)+omega*feqeq[k];
				
	}
	
}




void collide()
{
	for (int counter=0;counter<NUM;counter++)
		{
			rho[counter]=0.0;
		    
		    int offset=counter*NPOP;
	 
			double sum=0;
			for(int k=0;k<NPOP;k++)
				sum+=f[offset+k];
	
			rho[counter]=sum;

			double dense_temp=rho[counter];
			double ux_temp=ux[counter];
			double uy_temp=uy[counter];

			//BGK equilibrium
			double feq[NPOP];
	
			//TRT equilibrium
			double feq_plus[NPOP],feq_minus[NPOP];
		    double f_plus[NPOP],f_minus[NPOP];

		     //Speeding up the process
		    f_plus[0]=f[offset];
		    f_plus[1]=0.5*(f[offset+1]+f[offset+3]);
		    f_plus[2]=0.5*(f[offset+2]+f[offset+4]);
		    f_plus[3]=f_plus[1];
		    f_plus[4]=f_plus[2];
		    f_plus[5]=0.5*(f[offset+5]+f[offset+7]);
		    f_plus[6]=0.5*(f[offset+6]+f[offset+8]);
		    f_plus[7]=f_plus[5];
		    f_plus[8]=f_plus[6];

		    f_minus[0]=0.0;
		    f_minus[1]=0.5*(f[offset+1]-f[offset+3]);
		    f_minus[2]=0.5*(f[offset+2]-f[offset+4]);
		    f_minus[3]=-f_minus[1];
		    f_minus[4]=-f_minus[2];
		    f_minus[5]=0.5*(f[offset+5]-f[offset+7]);
		    f_minus[6]=0.5*(f[offset+6]-f[offset+8]);
		    f_minus[7]=-f_minus[5];
		    f_minus[8]=-f_minus[6];
		  
		    //The equilibrium populations are taken from Irina's bb_pressure and bb_examples
		    feq_minus[0]=0.0;
		    feq_minus[1]=weights_trt[1]*dense_temp*(cx[1]*ux_temp+cy[1]*uy_temp);
		    feq_minus[2]=weights_trt[2]*dense_temp*(cx[2]*ux_temp+cy[2]*uy_temp);
		    feq_minus[3]=-feq_minus[1];
		    feq_minus[4]=-feq_minus[2];
		    feq_minus[5]=weights_trt[5]*dense_temp*(cx[5]*ux_temp+cy[5]*uy_temp);
		    feq_minus[6]=weights_trt[6]*dense_temp*(cx[6]*ux_temp+cy[6]*uy_temp);
		    feq_minus[7]=-feq_minus[5];
		    feq_minus[8]=-feq_minus[6];
		    
		    double u_sq=ux_temp*ux_temp+uy_temp*uy_temp;
		
			//Check with Irina - need to correct t_q^{(u)}\neq t_q^{(a)}        
		    feq_plus[1]=weights_trt[1]*dense_temp*(ce+0.5*(3.0*(cx[1]*ux_temp+cy[1]*uy_temp)*(cx[1]*ux_temp+cy[1]*uy_temp)-u_sq));
		    feq_plus[2]=weights_trt[2]*dense_temp*(ce+0.5*(3.0*(cx[2]*ux_temp+cy[2]*uy_temp)*(cx[2]*ux_temp+cy[2]*uy_temp)-u_sq));
		    feq_plus[3]=feq_plus[1];
		    feq_plus[4]=feq_plus[2];
		    feq_plus[5]=weights_trt[5]*dense_temp*(ce+0.5*(3.0*(cx[5]*ux_temp+cy[5]*uy_temp)*(cx[5]*ux_temp+cy[5]*uy_temp)-u_sq));
		    feq_plus[6]=weights_trt[6]*dense_temp*(ce+0.5*(3.0*(cx[6]*ux_temp+cy[6]*uy_temp)*(cx[6]*ux_temp+cy[6]*uy_temp)-u_sq));
		    feq_plus[7]=feq_plus[5];
		    feq_plus[8]=feq_plus[6];
		    feq_plus[0]=dense_temp-2.0*(feq_plus[1]+feq_plus[2]+feq_plus[5]+feq_plus[6]);

	        if (dense_temp>0.5) 
                omega=omega_low;
            else
                omega=omega_high;
            omega_minus=omega;
            omega_plus=2.0-omega_minus;
			//Collision operator
			for(int k=0; k < NPOP; k++)
				f2[offset+k]=f[offset+k]-omega_plus*(f_plus[k]-feq_plus[k])-omega_minus*(f_minus[k]-feq_minus[k]);
	}
			
}



void initialize_geometry()
{
	NY=129;
	NX=129;
	NUM=NX*NY;
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    radius=40;
    
    for(int iX=0;iX<NX;iX++)
    	for(int iY=0;iY<NY;iY++)
    	{
    		int counter=iX*NY+iY;
    		if ((iX-(NX-1)/2)*(iX-(NX-1)/2)+(iY-(NY-1)/2)*(iY-(NY-1)/2)<=radius*radius)
    		{
    	    	rho[counter]=1;
                ux[counter]=0.0;
                uy[counter]=0.0;
    	    }
    	    else 
    	    {
    	    	rho[counter]=0.0;
                ux[counter]=0.0;
                uy[counter]=0.0;
    	    }
    	}
	
}


void finish_simulation()
{
	delete[] rho;
	delete[] ux;
	delete[] uy;
	delete[] f;
	delete[] f2;
}


void stream()
{
	for(int counter=0;counter<NUM;counter++)
	{	
        int iX=counter%NX;
        int iY=counter/NX;
        
		for(int iPop=0;iPop<NPOP;iPop++)
		{
            int iX2=(iX-cx[iPop]+NX)%NX;
            int iY2=(iY-cy[iPop]+NY)%NY;
			int counter2=iY2*NX+iX2;
			f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
		}
    }
}


int main(int argc, char* argv[])
{


    initialize_geometry();
    init();
    

	for(int counter=0;counter<=N;counter++)
	{

        collide();
		stream();
        
		//Writing files
		if (counter%NOUTPUT==0)
		{
			std::cout<<"Counter="<<counter<<"\n";
  			std::stringstream filewritedensity;
 			
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;

			filewritedensity<<"density"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			
 			writedensity(filewritedensity.str());
		}

	}

    finish_simulation();
   
   	return 0;
}
