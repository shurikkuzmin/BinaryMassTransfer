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

//Other constants
const int NPOP=9;

//Time steps
int N=50000;
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
double conc_inlet=0.5;
double u0=0.05;
std::vector<int> bb_nodes_inlet;
std::vector<char>* dirs_inlet;

std::vector<int> bb_nodes_wall;
std::vector<int> dirs_wall;

//BGK relaxation parameter
double omega=1.4; ///1.99;
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





void collide()
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

	 
			//Collision operator
			for(int k=0; k < NPOP; k++)
				f2[offset+k]=f[offset+k]-omega_plus*(f_plus[k]-feq_plus[k])-omega_minus*(f_minus[k]-feq_minus[k]);
	}
			
}


void update_bounce_back()
{
	//Updating the inlet populations through the anti bounce back
	for(int iX=2;iX<NX-2;iX++)
	{
		
		f2[iX*NPOP+2]=-f2[(NX+iX)*NPOP+4]+2*weights[2]*conc_inlet*(1.0+3.0*(ux[iX]*cx[2]+uy[iX]*cy[2]));
		f2[iX*NPOP+5]=-f2[(NX+iX+1)*NPOP+7]+2*weights[5]*conc_inlet*(1.0+3.0*(ux[iX]*cx[5]+uy[iX]*cy[5]));
		f2[iX*NPOP+6]=-f2[(NX+iX-1)*NPOP+8]+2*weights[6]*conc_inlet*(1.0+3.0*(ux[iX]*cx[6]+uy[iX]*cy[6]));
	}
	
	//Corners
	f2[NPOP+2]=-f2[(NX+1)*NPOP+4]+2*weights[2]*conc_inlet*(1.0+3.0*(ux[1]*cx[2]+uy[1]*cy[2]));
	f2[NPOP+5]=-f2[(NX+2)*NPOP+7]+2*weights[5]*conc_inlet*(1.0+3.0*(ux[1]*cx[5]+uy[1]*cy[5]));
	f2[(NX-2)*NPOP+2]=-f2[(2*NX-2)*NPOP+4]+2*weights[2]*conc_inlet*(1.0+3.0*(ux[NX-2]*cx[2]+uy[NX-2]*cy[2]));
	f2[(NX-2)*NPOP+6]=-f2[(2*NX-2)*NPOP+7]+2*weights[5]*conc_inlet*(1.0+3.0*(ux[NX-2]*cx[5]+uy[NX-2]*cy[5]));

    //Updating the outlet symmetry 
    for(int iX=1;iX<NX-1;iX++)
    {
    	for(int iPop=0;iPop<NPOP;iPop++)
    		f2[((NY-1)*NX+iX)*NPOP+iPop]=f2[((NY-2)*NX+iX)*NPOP+iPop];
    	rho[(NY-1)*NX+iX]=rho[(NY-2)*NX+iX];
    } 	

	
	
	//Updating the wall
	for(int iY=1;iY<NY-1;iY++)
	{
		f2[iY*NX*NPOP+1]=-f2[(iY*NX+1)*NPOP+3]+2*weights[1]*conc_wall*(1.0+3.0*(ux[iY*NX]*cx[1]+uy[iY*NX]*cy[1]));
		f2[iY*NX*NPOP+8]=-f2[((iY-1)*NX+1)*NPOP+6]+2*weights[8]*conc_wall*(1.0+3.0*(ux[iY*NX]*cx[8]+uy[iY*NX]*cy[8]));
		f2[iY*NX*NPOP+5]=-f2[((iY+1)*NX+1)*NPOP+7]+2*weights[5]*conc_wall*(1.0+3.0*(ux[iY*NX]*cx[5]+uy[iY*NX]*cy[5]));
	}
	
	//Corners
	f2[1]=-f2[NPOP+3]+2*weights[1]*conc_wall*(1.0+3.0*(ux[0]*cx[1]+uy[0]*cy[1]));
	f2[5]=-f2[(NX+1)*NPOP+7]+2*weights[5]*conc_wall*(1.0+3.0*(ux[0]*cx[5]+uy[0]*cy[5]));
	f2[(NY-1)*NX*NPOP+1]=-f2[((NY-1)*NX+1)*NPOP+3]+2*weights[1]*conc_wall*(1.0+3.0*(ux[(NY-1)*NX]*cx[1]+uy[(NY-1)*NX]*cy[1]));
	f2[(NY-1)*NX*NPOP+7]=-f2[((NY-2)*NX+1)*NPOP+5]+2*weights[7]*conc_wall*(1.0+3.0*(ux[(NY-1)*NX]*cx[7]+uy[(NY-1)*NX]*cy[7]));

	
	//Updating top wall symmetry
    for(int iY=0;iY<NY;iY++)
    {
    	for(int iPop=0;iPop<NPOP;iPop++)
    		f2[(iY*NX+NX-1)*NPOP+iPop]=f2[(iY*NX+NX-2)*NPOP+iPop];
    	rho[iY*NX+NX-1]=rho[iY*NX+NX-2];
    } 	
		
}

void initialize_geometry()
{
	NY=10*40;
	NX=40;
	NUM=NX*NY;
    geometry=new int[NUM];
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    
 
	//Initialization
    for(int iY=0;iY<NY;iY++)
    	for(int iX=0;iX<NX-1;iX++)
    	{
    		int counter=iY*NX+iX;
    		rho[counter]=0.0;
    		uy[counter]=u0*(1.0-double(iX*iX)/((NX-1.5)*(NX-1.5)));
    		ux[counter]=0.0;
    	}
 
    for(int iX=1;iX<NX-1;iX++)
    { 
    	rho[iX]=conc_inlet;
    	rho[(NY-1)*NX+iX]=rho[(NY-2)*NX+iX];
    }
    
    for(int iY=0;iY<NY;iY++)
    {
    	rho[iY*NX]=conc_wall;
        rho[iY*NX+NX-1]=rho[iY*NX+NX-2];
        ux[iY*NX+NX-1]=ux[iY*NX+NX-2];
        uy[iY*NX+NX-1]=ux[iY*NX+NX-2];
    }
	writedensity("conc_initial");
	writevelocityx("ux_initial");
	writevelocityy("uy_initial");    
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
				int counter2=counter-cy[iPop]*NX-cx[iPop];
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

        //collide(); 
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

			filewritedensity<<"film"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			
 			writedensity(filewritedensity.str());
		}

	}

    finish_simulation();
   
   	return 0;
}
