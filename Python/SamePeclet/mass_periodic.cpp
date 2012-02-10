//Simple mass transfer


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <limits>

//Domain size
int NY;
int NX;
int NUM;

//Overall scaling parameter
double peclet=1000.0;
double scale_velocity=1.0;
double scale_diffusion=1.0;

//Other constants
const int NPOP=9;

//Time steps
int N=1000000;
int NOUTPUT=20000;
int NSIGNAL=1000;

//Fields and populations
double *f;
double *f2;
double *rho;
double *ux;
double *uy;
int * geometry;

//Boundary conditions
double conc_inlet=0.0;
double conc_bubble=1.0;
double conc_wall=0.0;
std::vector<int> bb_nodes;
std::vector<char>* dirs;
int *bottom;
int *bottom_mid;
int *top;
int *top_mid;

//BGK relaxation parameter
double omega=1.99;
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

	for(int iY=0;iY<NY;iY++)
		for(int iX=0;iX<NX;iX++)
		{
			int  counter=iY*NX+iX;
			double dense_temp=rho[counter];
			double ux_temp=ux[counter];
			double uy_temp=uy[counter];
            
            sum=0;
            for (int k=0; k<NPOP; k++)
				sum+=f[counter*NPOP+k];
			
			rho[counter]=sum;
		}

	
	
	writedensity("check_initial_density");

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



void collide_column(int coor_y,int coor_bottom,int coor_top)
{
	for(int iX=coor_bottom;iX<=coor_top;iX++)
	{
		int counter=coor_y*NX+iX;
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
		double u_mn=ux_temp*ux_temp-uy_temp*uy_temp;
		double u_cr=ux_temp*uy_temp;
		
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

void collide_bulk()
{

    for(int iY=0;iY<NY;iY++)
		if (bottom[iY]==top[iY])
			collide_column(iY,1,NX-2);
			//collide_column_bgk(iY,1,NX-2);
		else
		{
			if(bottom_mid[iY]==top_mid[iY])
			{
				collide_column(iY,1,bottom[iY]);
				collide_column(iY,top[iY],NX-2);
			}
			else
			{
				collide_column(iY,1,bottom[iY]);
				collide_column(iY,bottom_mid[iY],top_mid[iY]);
				collide_column(iY,top[iY],NX-2);
			}
			
			//collide_column_bgk(iY,1,bottom[iY]);
			//collide_column_bgk(iY,top[iY],NX-2);
		}
			
}

void update_bounce_back()
{
	for(int counter=0;counter<bb_nodes.size();counter++)
	{
		for(int k=0;k<dirs[counter].size();k++)
		{
			int dir=dirs[counter][k];
			int counter2=bb_nodes[counter]+cy[dir]*NX+cx[dir];
			f2[bb_nodes[counter]*NPOP+dir]=-f2[counter2*NPOP+compliment[dir]]+2*weights[dir]*conc_bubble;
		}
		
	}
	
	//BB nodes density and velocity specification
	for(int iY=0;iY<NY;iY++)
	{
		int counter=iY*NX;
		int counter_top=((iY+1+NY)%NY)*NX;
		int counter_bottom=((iY-1+NY)%NY)*NX;


		f2[counter*NPOP+1]=f2[(counter+1)*NPOP+3];
		f2[counter*NPOP+5]=f2[(counter_top+1)*NPOP+7];
		f2[counter*NPOP+8]=f2[(counter_bottom+1)*NPOP+6];

		f2[(counter+NX-1)*NPOP+3]=f2[(counter+NX-2)*NPOP+1];
		f2[(counter+NX-1)*NPOP+6]=f2[(counter_top+NX-2)*NPOP+8];
		f2[(counter+NX-1)*NPOP+7]=f2[(counter_bottom+NX-2)*NPOP+5];

		rho[counter]=conc_wall;
		rho[counter+NX-1]=conc_wall;
		ux[counter]=0.0;
		ux[counter+NX-1]=0.0;
		uy[counter]=0.0;
		uy[counter+NX-1]=0.0;
	}

}

void initialize_geometry()
{
	NY=3000;
	NX=202;
	NUM=NX*NY;
    geometry=new int[NUM];
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    bottom=new int[NY];
    bottom_mid=new int[NY];
    top=new int[NY];
    top_mid=new int[NY];
    
	std::ifstream fin("geometry.dat");
	std::ifstream fux("vely0200000.dat");
	std::ifstream fuy("velx0200000.dat");
	
	//Reading files
	for(int counter=0;counter<NUM;counter++)
	{
		fin>>geometry[counter];
		fux>>uy[counter];
		fuy>>ux[counter];
	}
 
   	scale_velocity=0.05/uy[0];
  	scale_diffusion=0.05*200.0/peclet;
    omega=1.0/(3.0*scale_diffusion+0.5);
    omega_minus=omega;
    omega_plus=2.0-omega;
  
    

    //Scaling of velocities
    for(int counter=0;counter<NUM;counter++)
    {
    	ux[counter]=scale_velocity*ux[counter];
    	uy[counter]=scale_velocity*uy[counter];
    }


	//Initialization
    for(int counter=0;counter<NUM;counter++)
    {
		if (geometry[counter]==0)
		{
		    rho[counter]=conc_bubble;
			ux[counter]=0.0;
			uy[counter]=0.0;
			bb_nodes.push_back(counter);
		}
		else if(geometry[counter]==-1)
		{
			rho[counter]=-1.0;
			ux[counter]=0.0;
			uy[counter]=0.0;
		}
		else
			rho[counter]=0.0;
	}
	
	//Finding bulk nodes
	int iX=0;
	int iY=0;
	while(iY<NY)
	{
	    bool flag=false;
	    iX=0;
	    while(iX<NX)
	    {
			int counter=iY*NX+iX;
		    if (geometry[counter]==0) 
		    {
				flag=true;
				bottom[iY]=iX-1;
				while((geometry[counter]==0) || (geometry[counter]==-1))
				{
					iX++;
					counter=iY*NX+iX;
				}
				if (iX<=NX/2)
				{
					bottom_mid[iY]=iX;
					while(geometry[counter]==1)
					{
						iX++;
						counter=iY*NX+iX;
					}
					top_mid[iY]=iX-1;
					while((geometry[counter]==0) || (geometry[counter]==-1))
					{
						iX++;
						counter=iY*NX+iX;
					}
				    top[iY]=iX;
	
				}
				else
				{
					top[iY]=iX;
					bottom_mid[iY]=top_mid[iY]=0;
				}
			}
			iX++;
		}
		if (!flag)
			bottom[iY]=bottom_mid[iY]=top[iY]=top_mid[iY]=0;
		iY++;
	}
	
	std::ofstream fcoor("coor.dat");
	for(int iY=0;iY<NY;iY++)
	{
		fcoor<<bottom[iY]<<" ";
	}
	fcoor<<"\n";
	for(int iY=0;iY<NY;iY++)
	{
		fcoor<<top[iY]<<" ";
	}
	fcoor<<"\n";
	for(int iY=0;iY<NY;iY++)
	{
		fcoor<<bottom_mid[iY]<<" ";
	}
	fcoor<<"\n";
	for(int iY=0;iY<NY;iY++)
	{
		fcoor<<top_mid[iY]<<" ";
	}
	
	//Finding directions for BB nodes
    dirs=new std::vector<char>[bb_nodes.size()];
    for(int counter=0;counter<bb_nodes.size();counter++)
	{
		for(int k=1;k<NPOP;k++)
		{
			int counter2=bb_nodes[counter]+cy[k]*NX+cx[k];
			if (geometry[counter2]==1)
				dirs[counter].push_back(k);
		}
	}
}


void finish_simulation()
{
	delete[] geometry;
	delete[] rho;
	delete[] ux;
	delete[] uy;
	delete[] f;
	delete[] f2;
	delete[] bottom;
	delete[] bottom_mid;
	delete[] top_mid;
	delete[] top;
	delete[] dirs;
}

void stream_column(int coor_y,int coor_bottom,int coor_top)
{

	for(int iX=coor_bottom;iX<=coor_top;iX++)
	{
		int counter=coor_y*NX+iX;
		for(int iPop=0;iPop<NPOP;iPop++)
		{
			int iX2=(iX-cx[iPop]+NX)%NX;
			int iY2=(coor_y-cy[iPop]+NY)%NY;
						
			int counter2=iY2*NX+iX2;
			f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
		}
	}
}

void stream()
{
	
    for(int iY=0;iY<NY;iY++)
		if (bottom[iY]==top[iY])
			stream_column(iY,1,NX-2);
		else
		{
			if(top_mid[iY]==bottom_mid[iY])
			{
				stream_column(iY,1,bottom[iY]);
				stream_column(iY,top[iY],NX-2);
			}
			else
			{
				stream_column(iY,1,bottom[iY]);
				stream_column(iY,bottom_mid[iY],top_mid[iY]);
				stream_column(iY,top[iY],NX-2);
			}
		}

	
}

void calculate_mass_transfer(int time_counter)
{
	static std::ofstream fconc("concentration_periodic.dat");
	
    double conc=0.0;
    int counter; 	
    for(int iY=0;iY<NY;iY++)
		if (bottom[iY]==top[iY])
			for(int iX=1;iX<NX-1;iX++)
				conc+=rho[iY*NX+iX];
		else
		{
			if (top_mid[iY]==bottom_mid[iY])
			{
				for(int iX=1;iX<=bottom[iY];iX++)
					conc+=rho[iY*NX+iX];
				for(int iX=top[iY];iX<NX-1;iX++)
					conc+=rho[iY*NX+iX];
			}
			else
			{
				for(int iX=1;iX<=bottom[iY];iX++)
					conc+=rho[iY*NX+iX];
				for(int iX=bottom_mid[iY];iX<=top_mid[iY];iX++)
					conc+=rho[iY*NX+iX];
				for(int iX=top[iY];iX<NX-1;iX++)
					conc+=rho[iY*NX+iX];
			}
		}
		
	double conc_outlet=0.0;
	double sum_vel=0.0;
	for(int iX=1;iX<NX-1;iX++)
	{
		conc_outlet+=rho[iX]*uy[iX];
		sum_vel+=uy[iX];
	}	

	fconc<<time_counter<<" "<<conc<<" "<<conc_outlet/sum_vel<<"\n"<<std::flush;
}

int main(int argc, char* argv[])
{

    if (argc==2)
    	peclet=atof(argv[1]);
    	
    initialize_geometry();
    init();
 
    std::cout<<"Scaling velocity parameter is "<<scale_velocity<<"\n";
    std::cout<<"Omega="<<omega<<"\n";
    std::cout<<"Peclet="<<peclet<<"\n";
      

	for(int counter=0;counter<=N;counter++)
	{

        collide_bulk();
        update_bounce_back();
		stream();
        
        if(counter%NSIGNAL==0)
        {
        	if (rho[NX/2]!=rho[NX/2])
			{
				std::cout<<"Simulation is unstable!\n";
				return 0;
			}
		    calculate_mass_transfer(counter);
       
        }
		//Writing files
		if (counter%NOUTPUT==0)
		{
			std::cout<<"Counter="<<counter<<"\n";
  			std::stringstream filewritedensity;
 			
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;

			filewritedensity<<"density_periodic"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			
 			writedensity(filewritedensity.str());
		}

	}

    finish_simulation();
   
   	return 0;
}
