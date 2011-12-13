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

//Other constants
const int NPOP=9;

//Time steps
int N=50000;
int NOUTPUT=1000;

//Fields and populations
double *f;
double *f2;
double *rho;
double *ux;
double *uy;
int * geometry;

//Boundary conditions
double conc_wall=1.0;
double conc_inlet=0.0;
double u0=0.025;
std::vector<int> bb_nodes_inlet;
std::vector<char>* dirs_inlet;

std::vector<int> bb_nodes_wall;
std::vector<int> dirs_wall;

//BGK relaxation parameter
double omega=1.8; ///1.99;
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

void writevelocityx(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iX=0; iX<NX; ++iX)
	{
		for (int iY=0; iY<NY; iY++)
		{
			int counter=iY*NX+iX;
			fout<<ux[counter]<<" ";
		}
		fout<<"\n";
	}
}

void writevelocityy(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iX=0; iX<NX; ++iX)
	{
		for (int iY=0; iY<NY; iY++)
		{
			int counter=iY*NX+iX;
			fout<<uy[counter]<<" ";
		}
		fout<<"\n";
	}
}


void initialize_geometry()
{
	NY=20*80;
	NX=80;
	NUM=NX*NY;
    geometry=new int[NUM];
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    
 
	//Initialization
    for(int iY=0;iY<NY;iY++)
    	for(int iX=0;iX<NX;iX++)
    	{
    		int counter=iY*NX+iX;
    		rho[counter]=0.0;
    		ux[counter]=0.0;
    		double factor=double(iX-0.5)/double(NX-1.5);
    		uy[counter]=4.0*u0*(1.0-factor)*factor;
    		if ((iX==0) || (iX==NX-1))
    			uy[counter]=0.0;
    	}
 
    for(int iY=0;iY<NY;iY++)
    {
    	rho[iY*NX]=conc_wall;
        rho[iY*NX+NX-1]=conc_wall;
    }
	writedensity("conc_initial");
	writevelocityx("ux_initial");
	writevelocityy("uy_initial");    
}



void init()
{
	//Creating arrays
	f=new double[NUM*NPOP];
	f2=new double[NUM*NPOP];
	
	//Bulk nodes initialization
	double feq;
	
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


void collide_bgk()
{
	for (int iY=1;iY<NY-1;iY++)
		for(int iX=1;iX<NX-1;iX++)
		{
			int counter=iY*NX+iX;
			rho[counter]=0.0;
		    
		    int offset=counter*NPOP;
	 
			double sum=0;
			for(int k=0;k<NPOP;k++)
				sum+=f[offset+k];
	
			rho[counter]=sum;

			double dense_temp=rho[counter];
			double ux_temp=ux[counter];
			double uy_temp=uy[counter];
			
			double feq[NPOP];
	 
			//Collision operator
			for(int k=0; k < NPOP; k++)
			{
				feq[k]=weights[k]*dense_temp*(1.0+3.0*(cx[k]*ux_temp+cy[k]*uy_temp)+4.5*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+2.0*cx[k]*cy[k]*ux_temp*uy_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp));
				f2[offset+k]=f[offset+k]-omega*(f[offset+k]-feq[k]);
			}
	}
			
}


void update_bounce_back()
{

	//Updating outlet symmetry
    for(int iX=1;iX<NX-1;iX++)
    {
    	for(int iPop=0;iPop<NPOP;iPop++)
		{
    		f2[((NY-1)*NX+iX)*NPOP+iPop]=f2[((NY-2)*NX+iX)*NPOP+iPop];
    	}
    	rho[(NY-1)*NX+iX]=rho[(NY-2)*NX+iX];
    } 	
	
	//Updating walls
	for(int iY=1;iY<NY-1;iY++)
	{
		f2[iY*NX*NPOP+1]=-f2[(iY*NX+1)*NPOP+3]+2.0*weights[1]*conc_wall;
		f2[iY*NX*NPOP+8]=-f2[((iY-1)*NX+1)*NPOP+6]+2.0*weights[8]*conc_wall;
		f2[iY*NX*NPOP+5]=-f2[((iY+1)*NX+1)*NPOP+7]+2.0*weights[5]*conc_wall;

		f2[(iY*NX+NX-1)*NPOP+3]=-f2[(iY*NX+NX-2)*NPOP+1]+2.0*weights[3]*conc_wall;
		f2[(iY*NX+NX-1)*NPOP+6]=-f2[((iY+1)*NX+NX-2)*NPOP+8]+2.0*weights[6]*conc_wall;
		f2[(iY*NX+NX-1)*NPOP+7]=-f2[((iY-1)*NX+NX-2)*NPOP+5]+2.0*weights[7]*conc_wall;
		rho[iY*NX]=conc_wall;
		rho[iY*NX+NX-1]=conc_wall;
	}
	
	//Updating inlet
	for(int iX=1;iX<NX-1;iX++)
	{
		f2[iX*NPOP+2]=-f2[(NX+iX)*NPOP+4]+2.0*weights[2]*conc_inlet;
		f2[iX*NPOP+6]=-f2[(NX+iX-1)*NPOP+8]+2.0*weights[6]*conc_inlet;
		f2[iX*NPOP+5]=-f2[(NX+iX+1)*NPOP+7]+2.0*weights[5]*conc_inlet;
		rho[iX]=conc_inlet;
	}
	

	//Corners
	f2[1]=-f2[NPOP+3]+2.0*weights[1]*conc_wall;
	f2[5]=-f2[(NX+1)*NPOP+7]+2.0*weights[5]*conc_wall;
	f2[(NY-1)*NX*NPOP+1]=-f2[((NY-1)*NX+1)*NPOP+3]+2.0*weights[1]*conc_wall;
	f2[(NY-1)*NX*NPOP+8]=-f2[((NY-2)*NX+1)*NPOP+6]+2.0*weights[7]*conc_wall;

	f2[(NX-1)*NPOP+3]=-f2[(NX-2)*NPOP+1]+2.0*weights[3]*conc_wall;
	f2[(NX-1)*NPOP+6]=-f2[(NX+NX-2)*NPOP+8]+2.0*weights[6]*conc_wall;
	f2[((NY-1)*NX+NX-1)*NPOP+3]=-f2[((NY-1)*NX+NX-2)*NPOP+1]+2.0*weights[3]*conc_wall;
	f2[((NY-1)*NX+NX-1)*NPOP+7]=-f2[((NY-2)*NX+NX-2)*NPOP+5]+2.0*weights[7]*conc_wall;


	rho[0]=conc_wall;
	rho[NX-1]=conc_wall;
	rho[(NY-1)*NX]=conc_wall;
	rho[(NY-1)*NX+NX-1]=conc_wall;
		
}



void finish_simulation()
{
	delete[] geometry;
	delete[] rho;
	delete[] ux;
	delete[] uy;
	delete[] f;
	delete[] f2;
}


void stream()
{
	for (int iY=1;iY<NY-1;iY++)
		for(int iX=1;iX<NX-1;iX++)
		{
			int counter=iY*NX+iX;
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

        collide_bgk();
        update_bounce_back();
		stream();
        
		//Writing files
		if (counter%NOUTPUT==0)
		{
			std::cout<<"Counter="<<counter<<"\n";
  			std::stringstream filewritedensity;
 			
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;

			filewritedensity<<"film_antibb"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			
 			writedensity(filewritedensity.str());
		}

	}

    finish_simulation();
   
   	return 0;
}
