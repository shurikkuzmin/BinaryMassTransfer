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
	std::ifstream fux("ux_non.dat");
	std::ifstream fuy("uy_non.dat");
	std::ifstream fdiv("divergence.dat");
	
	//Reading files
	for(int counter=0;counter<NUM;counter++)
	{
		fin>>geometry[counter];
		fux>>uy[counter];
		fuy>>ux[counter];
		fdiv>>divergence[counter];
	}
