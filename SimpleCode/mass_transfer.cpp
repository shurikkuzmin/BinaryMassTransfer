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
int N=40000;
int NOUTPUT=1000;

//Fields and populations
double *f;
double *f2;
double *rho;
double *ux;
double *uy;
int * geometry;

//double f[NX][NY][9], f2[NX][NY][9], g[NX][NY][9], g2[NX][NY][9];
//double rho[NX][NY],ux[NX][NY],uy[NX][NY],phase[NX][NY];

//Boundary conditions
double conc_inlet=0.5;
double conc_bubble=1.0;

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
	std::ofstream fout(fname.c_str());
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


//void collide_bulk()
//{
    ////The phase field should be calculated prior the laplacians
    //for(int iX=0;iX<NX;iX++)
        //for(int iY=1;iY<NY-1;iY++)
		//{
            //phase[iX][iY]=0.0;
            //for(int iPop=0;iPop<9;iPop++)
   				//phase[iX][iY]+=g[iX][iY][iPop];
		//}

    ////The phase of the BB nodes
    //for(int iX=0;iX<NX;iX++)
    //{
		////First-order accuracy
		//phase[iX][0]=phase[iX][1]-wall_gradient;
        ////Second-order accuracy
        ////phase[iX][iY]=(4.0*phase[iX][iY+1]-phase[iX][iY+2]-2.0*wall_gradient)/3.0;
        ////First-order accuracy
		//phase[iX][NY-1]=phase[iX][NY-2]-wall_gradient;
		////Second-order accuracy
        ////phase[iX][iY]=(4.0*phase[iX][iY-1]-phase[iX][iY-2]-2.0*wall_gradient)/3.0;
    //}

    //for(int iX=0;iX<NX;iX++)
        //for(int iY=1;iY<NY-1;iY++)
		//{

			////Construction equilibrium
			//rho[iX][iY]=0.0;
			//ux[iX][iY]=0.0;
			//uy[iX][iY]=0.0;

			//for(int iPop=0;iPop<9;iPop++)
			//{
				//rho[iX][iY]+=f[iX][iY][iPop];
				//ux[iX][iY]+=f[iX][iY][iPop]*cx[iPop];
				//uy[iX][iY]+=f[iX][iY][iPop]*cy[iPop];
			//}

			//ux[iX][iY]=(ux[iX][iY]+0.5*force_x)/rho[iX][iY];
			//uy[iX][iY]=(uy[iX][iY]+0.5*force_y)/rho[iX][iY];


			//double laplace_temp=0.0;
			//double gradx_temp=0.0;
			//double grady_temp=0.0;
			//for(int k=0;k<9;k++)
			//{
				//int iX2=(iX+cx[k]+NX) % NX;
				//int iY2=(iY+cy[k]+NY) % NY;
				//laplace_temp+=laplacestencil[k]*phase[iX2][iY2];
				//gradx_temp+=gradstencilx[k]*phase[iX2][iY2];
				//grady_temp+=gradstencily[k]*phase[iX2][iY2];
			//}

			//double phase_temp=phase[iX][iY];
			//double dense_temp=rho[iX][iY];
			//double ux_temp=ux[iX][iY];
			//double uy_temp=uy[iX][iY];

			//double sum=0.0;
			//double sum_phase=0.0;
			//double phase_square=phase_temp*phase_temp;
			//double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
			//double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

			//double feqeq[9],geqeq[9],force[9];

			////Obtain force population
            //for (int k=0;k<9;k++)
            //{
				//force[k]=weights[k]*(1.0-0.5*omega)*(3.0*force_x*cx[k]+3.0*force_y*cy[k]+
                            //9.0*((cx[k]*cx[k]-1.0/3.0)*force_x*ux_temp+cx[k]*cy[k]*(force_x*uy_temp+force_y*ux_temp)+
 							//(cy[k]*cy[k]-1.0/3.0)*force_y*uy_temp));
            //}
			//for (int k=1; k<9; k++)
			//{
				//feqeq[k]=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								//+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
				//+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
				//geqeq[k]=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								//+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				//sum+=feqeq[k];
				//sum_phase+=geqeq[k];

			//}

			//feqeq[0]=dense_temp-sum;
			//geqeq[0]=phase_temp-sum_phase;

            //double tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
            //double omega_temp=1.0/tau_temp;


			//for(int k=0; k < 9; k++)
			//{
				//f2[iX][iY][k]=f[iX][iY][k]*(1.0-omega_temp)+omega_temp*feqeq[k]+force[k];
				//g2[iX][iY][k]=g[iX][iY][k]*(1.0-omega_temp)+omega_temp*geqeq[k];
			//}


		//}

//}

//void update_bounce_back()
//{
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

//}

void initialize_geometry()
{
	NY=202;
	NX=3001;
	NUM=NX*NY;
    geometry=new int[NUM];
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
	std::ifstream fin("geometry.dat");
	std::ifstream fux("ux.dat");
	std::ifstream fuy("uy.dat");
	
	for(int counter=0;counter<NUM;counter++)
	{
		fin>>geometry[counter];
		fux>>ux[counter];
		fuy>>uy[counter];
	}

    for(int counter=0;counter<NUM;counter++)
    {
		if (geometry[counter]==0)
		{
		    rho[counter]=conc_bubble;
			ux[counter]=0;
			uy[counter]=0;
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
	//writedensity("conc.dat");
    
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

int main(int argc, char* argv[])
{

    initialize_geometry();
    init();
    finish_simulation();
    //if (argc!=1)
    //{
        //width=atoi(argv[1]);
	//int ratio=atoi(argv[2]);
	//force_x=0.000006/(ratio*ratio);
	////NY=49*ratio+2;
	////NX=49*ratio*25+1;
	//N=40000*ratio;
	//NOUTPUT=1000*ratio;
    //}
    //std::cout<<"Width="<<width<<"\n";

    //matrix_init();
    //init();


	//for(int counter=0;counter<=N;counter++)
	//{

        //collide_bulk();
        //update_bounce_back();

		////Streaming
		//for(int iX=0;iX<NX;iX++)
			//for(int iY=1;iY<NY-1;iY++)
				//for(int iPop=0;iPop<9;iPop++)
				//{
					//int iX2=(iX-cx[iPop]+NX)%NX;
					//int iY2=(iY-cy[iPop]+NY)%NY;
					//f[iX][iY][iPop]=f2[iX2][iY2][iPop];
                    //g[iX][iY][iPop]=g2[iX2][iY2][iPop];
				//}

		////Writing files
		//if (counter%NOUTPUT==0)
		//{
			//std::cout<<counter<<"\n";

			//std::stringstream filewritedensity;
 			//std::stringstream filewritevelocity;
 			//std::stringstream filewritephase;
 			//std::stringstream counterconvert;
 			//counterconvert<<counter;
 			//filewritedensity<<std::fixed;
			//filewritevelocity<<std::fixed;
			//filewritephase<<std::fixed;

			//filewritedensity<<"density"<<std::string(6-counterconvert.str().size(),'0')<<counter;
			//filewritevelocity<<"velocity"<<std::string(6-counterconvert.str().size(),'0')<<counter;
            //filewritephase<<"phase"<<std::string(6-counterconvert.str().size(),'0')<<counter;

 			//writedensity(filewritedensity.str());
			//writevelocity(filewritevelocity.str());
            //writephase(filewritephase.str());
		//}


	//}

	//return 0;
}
