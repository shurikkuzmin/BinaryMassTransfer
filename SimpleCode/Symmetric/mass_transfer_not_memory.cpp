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
int N=10000;
int NOUTPUT=1000;

//Fields and populations
double *f;
double *f2;
double *rho;
double *ux;
double *uy;
int * geometry;

//Boundary conditions
double conc_inlet=0.5;
double conc_bubble=1.0;
std::vector<int> bb_nodes;
std::vector<char>* dirs;
int *bottom;
int *top;

//BGK relaxation parameter
double omega=1.0;
double omega_plus=omega;
double omega_minus=omega;

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

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
		{
			int counter=iY*NX+iX;
			fout<<rho[counter]<<" ";
		}
		fout<<"\n";
	}

}

//void writevelocityx(std::string const & fName)
//{
	//std::string fullName = "./tmp2/" + fName+ ".dat";
	//std::ofstream fout(fullName.c_str());
	//fout.precision(10);

	//for (int iY=NY-1; iY>=0; --iY)
	//{
		//for (int iX=0; iX<NX; ++iX)
			//fout<<ux[iX][iY]<<" ";
		//fout<<"\n";
	//}

//}

//void writevelocityy(std::string const & fName)
//{
	//std::string fullName = "./tmp2/" + fName+ ".dat";
	//std::ofstream fout(fullName.c_str());
	//fout.precision(10);

	//for (int iY=NY-1; iY>=0; --iY)
	//{
		//for (int iX=0; iX<NX; ++iX)
			//fout<<uy[iX][iY]<<" ";
		//fout<<"\n";
	//}

//}

void init()
{
	//Creating arrays
	f=new double[NUM*NPOP];
	f2=new double[NUM*NPOP];
	
	//Bulk nodes initialization
	double feq;
	double geq;
	double sum;
	
	for(int iX=0;iX<NX;iX++)
		for(int iY=0;iY<NY;iY++)
		{
			int  counter=iY*NX+iX;
			double dense_temp=rho[counter];
			double ux_temp=ux[counter];
			double uy_temp=uy[counter];
            
            sum=0;
            for (int k=0; k<NPOP; k++)
			{
				feq=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
								                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
								                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				sum+=feq;
				f[counter*NPOP+k]=feq;
			}
			
			//f[counter*NPOP]=dense_temp-sum;
		}

	for(int iX=0;iX<NX;iX++)
		for(int iY=0;iY<NY;iY++)
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

	
	
	writedensity("check_density.dat");


	////BB nodes initialization
    //for(int iX=0;iX<NX;iX++)
	//{
		//rho[iX][0]=1.0;
		//ux[iX][0]=0.0;
		//uy[iX][0]=0.0;

		//double fluxx=3.0*ux[iX][0];
		//double fluxy=3.0*uy[iX][0];

		//double qxx=4.5*ux[iX][0]*ux[iX][0];
		//double qxy=9.0*ux[iX][0]*uy[iX][0];
		//double qyy=4.5*uy[iX][0]*uy[iX][0];


		//for(int iPop=0;iPop<9;iPop++)
		//{
			//f[iX][0][iPop]=weights[iPop]*rho[iX][0]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
							//qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
			//g[iX][0][iPop]=weights[iPop]*phase[iX][0]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
														   //qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));

		//}

		//rho[iX][NY-1]=1.0;
        //ux[iX][NY-1]=0.0;
        //uy[iX][NY-1]=0.0;

		//fluxx=3.0*ux[iX][NY-1];
		//fluxy=3.0*uy[iX][NY-1];

		//qxx=4.5*ux[iX][NY-1]*ux[iX][NY-1];
		//qxy=9.0*ux[iX][NY-1]*uy[iX][NY-1];
		//qyy=4.5*uy[iX][NY-1]*uy[iX][NY-1];


		//for(int iPop=0;iPop<9;iPop++)
		//{
			//f[iX][NY-1][iPop]=weights[iPop]*rho[iX][NY-1]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
							//qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
			//g[iX][NY-1][iPop]=weights[iPop]*phase[iX][NY-1]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
														   //qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));

		//}
	//}

}

void collide_column(int coor_x,int coor_bottom,int coor_top)
{
	for(int iY=coor_bottom;iY<=coor_top;iY++)
	{
		int counter=iY*NX+coor_x;
		rho[counter]=0.0;

		double sum=0;
		for(int k=0;k<9;k++)
			sum+=f[counter*NPOP+k];
	
		rho[counter]=sum;

		double dense_temp=rho[counter];
		double ux_temp=ux[counter];
		double uy_temp=uy[counter];

		double feqeq[9];

		for (int k=0; k<9; k++)
			feqeq[k]=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
            						         +4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
			        				                         +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
					    			                         +2.0*ux_temp*uy_temp*cx[k]*cy[k]));

		for(int k=0; k < 9; k++)
			f2[counter*NPOP+k]=f[counter*NPOP+k]*(1.0-omega)+omega*feqeq[k];
				
	}
	
}

void collide_bulk()
{

    for(int iX=0;iX<NX;iX++)
		if (bottom[iX]==top[iX])
			collide_column(iX,0,NY-1);
		else
		{
			collide_column(iX,0,bottom[iX]);
			collide_column(iX,top[iX],NY-1);
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
	
	
	////BB nodes density and velocity specification
	//for(int iX=0;iX<NX;iX++)
	//{
		//int iXtop=(iX+1+NX)%NX;
		//int iXbottom=(iX-1+NX)%NX;

		//f2[iX][0][2]=f2[iX][1][4];
		//f2[iX][0][5]=f2[iXtop][1][7];
		//f2[iX][0][6]=f2[iXbottom][1][8];

		//f2[iX][NY-1][4]=f2[iX][NY-2][2];
		//f2[iX][NY-1][7]=f2[iXbottom][NY-2][5];
		//f2[iX][NY-1][8]=f2[iXtop][NY-2][6];

		////BB for the scalar phase??? Ask Irina about it
		//g2[iX][0][2]=g2[iX][1][4];
		//g2[iX][0][5]=g2[iXtop][1][7];
		//g2[iX][0][6]=g2[iXbottom][1][8];

		//g2[iX][NY-1][4]=g2[iX][NY-2][2];
		//g2[iX][NY-1][7]=g2[iXbottom][NY-2][5];
		//g2[iX][NY-1][8]=g2[iXtop][NY-2][6];


		//rho[iX][0]=1.0;
		//rho[iX][NY-1]=1.0;
		////phase[iX][0]=0.5;
		////phase[iX][NY-1]=0.5;
		//ux[iX][0]=0.0;ux[iX][NY-1]=0.0;
		//uy[iX][0]=0.0;uy[iX][NY-1]=0.0;
	//}

}

void initialize_geometry()
{
	NY=202;
	NX=3001;
	NUM=NX*NY;
    geometry=new int[NUM];
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    bottom=new int[NX];
    top=new int[NX];
    
	std::ifstream fin("geometry.dat");
	std::ifstream fux("ux.dat");
	std::ifstream fuy("uy.dat");
	
	//Reading files
	for(int counter=0;counter<NUM;counter++)
	{
		fin>>geometry[counter];
		fux>>ux[counter];
		fuy>>uy[counter];
	}

	//Initialization
    for(int counter=0;counter<NUM;counter++)
    {
		if (geometry[counter]==0)
		{
		    rho[counter]=conc_bubble;
			ux[counter]=0;
			uy[counter]=0;
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
	while(iX<NX)
	{
	    bool flag=false;
	    iY=0;
	    while(iY<NY)
	    {
			int counter=iY*NX+iX;
		    if (geometry[counter]==0) 
		    {
				flag=true;
				bottom[iX]=iY-1;
				while((geometry[counter]==0) || (geometry[counter]==-1))
				{
					iY++;
					counter=iY*NX+iX;
				}
				top[iX]=iY;
			}
			iY++;
		}
		if (!flag)
			bottom[iX]=top[iX]=0;
		iX++;
	}
	
	std::ofstream fcoor("coor.dat");
	for(int iX=0;iX<NX;iX++)
	{
		fcoor<<bottom[iX]<<" ";
	}
	fcoor<<"\n";
	for(int iX=0;iX<NX;iX++)
	{
		fcoor<<top[iX]<<" ";
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
	delete[] top;
	delete[] dirs;
}

void stream_column(int coor_x,int coor_bottom,int coor_top)
{

	for(int iY=coor_bottom;iY<=coor_top;iY++)
	{
		int counter=iY*NX+coor_x;
		for(int iPop=0;iPop<NPOP;iPop++)
		{
			int iX2=(coor_x-cx[iPop]+NX)%NX;
			int iY2=(iY-cy[iPop]+NY)%NY;
						
			int counter2=iY2*NX+iX2;
			f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
		}
	}
}

void stream()
{
	
    for(int iX=0;iX<NX;iX++)
		if (bottom[iX]==top[iX])
			stream_column(iX,0,NY-1);
		else
		{
			stream_column(iX,0,bottom[iX]);
			stream_column(iX,top[iX],NY-1);
		}

	
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

			filewritedensity<<"density"<<std::string(6-counterconvert.str().size(),'0')<<counter;
			
 			writedensity(filewritedensity.str());
		}

	}

    finish_simulation();
   
   	return 0;
}
