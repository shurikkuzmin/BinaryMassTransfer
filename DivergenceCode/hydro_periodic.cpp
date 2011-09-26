//Simple mass transfer


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
int N=20000;
int NOUTPUT=100;

//Fields and populations
double *f;
double *f2;
double *rho;
double *ux;
double *uy;
int * geometry;
double * divergence;

//Boundary conditions
double conc_inlet=0.5;
double conc_bubble=1.0;
double conc_wall=1.0;
std::vector<int> bb_nodes;
std::vector<char>* dirs;
int *bottom;
int *bottom_mid;
int *top;
int *top_mid;

//BGK relaxation parameter
double omega=2.5;

//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
int cx[]={0,1,0,-1,0,1,-1,-1,1};
int cy[]={0,0,1,0,-1,1,1,-1,-1};
int compliment[]={0,3,4,1,2,7,8,5,6};

void writedensity(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
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

void collide_bulk()
{

    for(int iY=0;iY<NY;iY++)
		if (bottom[iY]==top[iY])
			collide_column_bgk(iY,1,NX-2);
		else
		{
			if(bottom_mid[iY]==top_mid[iY])
			{
				collide_column_bgk(iY,1,bottom[iY]);
				collide_column_bgk(iY,top[iY],NX-2);
			}
			else
			{
				collide_column_bgk(iY,1,bottom[iY]);
				collide_column_bgk(iY,bottom_mid[iY],top_mid[iY]);
				collide_column_bgk(iY,top[iY],NX-2);
			}
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
    divergence=new double[NUM];
    
    bottom=new int[NY];
    bottom_mid=new int[NY];
    top=new int[NY];
    top_mid=new int[NY];
   
	std::ifstream fin("geometry.dat");
	std::ifstream fux("ux.dat");
	std::ifstream fuy("uy.dat");
	std::ifstream fdiv("divergence.dat");
	
	//Reading files
	for(int counter=0;counter<NUM;counter++)
	{
		fin>>geometry[counter];
		fux>>uy[counter];
		fuy>>ux[counter];
		fdiv>>divergence[counter];
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
			rho[counter]=1.0;
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
	static std::ofstream fconc("concentration.dat");
	
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

    initialize_geometry();
    init();
      

	for(int counter=0;counter<=N;counter++)
	{

        collide_bulk();
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

			filewritedensity<<"density"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			
 			writedensity(filewritedensity.str());
		    calculate_mass_transfer(counter);
		}

	}

    finish_simulation();
   
   	return 0;
}
