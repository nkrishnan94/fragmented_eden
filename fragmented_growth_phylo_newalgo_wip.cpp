#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <ctime> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <tuple>



///-----------SIMULATION PARAMETERS--------
const int nDemes = 501; /// side length of square simulation lattice
int initRad = 50; // initial innoculation radius in demes
int patchS = 20;
int patchP = 5;
unsigned int seedID = 37; //rng sseed
unsigned int nGen = 50000;


int checkEmpty(long double popA[nDemes][nDemes], int x, int y, const int arrSize ){

	int arr[2] = {x, y}; // store current lattice coordinates
	int neighbDiff[4][2] = {{0,1},{0,-1},{1,0},{-1,0}}; //directions for negihbors
	int neighbCoor[4][2]; // initialize array for neighbor coordinates
	std::vector <int> neighbEmpty;
	//int fullCount =0;
	int ne =0;
	int emptyFlag =0;


	while ((ne <4 ) && (emptyFlag ==0 )){

		neighbCoor[ne][0] = (arr[0] + nDemes+neighbDiff[ne][0]) % nDemes; 
		neighbCoor[ne][1] = (arr[1] + nDemes+neighbDiff[ne][1]) % nDemes;

		if (popA[ neighbCoor[ne][0] ][ neighbCoor[ne][1] ] == 0){

			emptyFlag=1;
		}
		ne+=1;
	}

	return emptyFlag; //returns 1 if an empty neighbor is present, otherwise 0;

}



std::vector <int> pickEmpty(long double popA[nDemes][nDemes], int x, int y, const int arrSize, const gsl_rng *R){


	int arr[2] = {x, y}; // store current lattice coordinates
	int neighbDiff[4][2] = {{0,1},{0,-1},{1,0},{-1,0}}; //directions for negihbors
	int neighbCoor[4][2]; // initialize array for neighbor coordinates
	std::vector <int> neighbEmpty;
	//int fullCount =0;





	for (int ne =0; ne< 4;ne ++){

		neighbCoor[ne][0] = (arr[0] + nDemes+neighbDiff[ne][0]) % nDemes; 
		neighbCoor[ne][1] = (arr[1] + nDemes+neighbDiff[ne][1]) % nDemes;
		

		if (popA[ neighbCoor[ne][0] ][ neighbCoor[ne][1] ] ==0 ){

			neighbEmpty.push_back(ne);


		}

	}



	if (neighbEmpty.size()==0 ){
		std::cout << "Deme supposed to have empty neighbor " << std::endl;
        exit(EXIT_FAILURE);
	}

	int nePick = gsl_rng_uniform_int(R, neighbEmpty.size() );


	std::vector <int> neighbCoorVec;
	neighbCoorVec.push_back(neighbCoor[ neighbEmpty[nePick] ] [0]);
	neighbCoorVec.push_back(neighbCoor[ neighbEmpty[nePick] ] [1] );

	


	return neighbCoorVec;

}



int main (int argc, char * argv[]){
	using namespace std;

	////-----DEFINE INPUT FLAGS----////

	int c;
    while ((c = getopt (argc, argv, "I:S:P:R")) != -1)
    {
        if (c == 'I')
            seedID  = atoi(optarg); // carrying capacity
        else if (c=='S')
        	patchS = atoi(optarg);
        else if (c=='P')
        	patchP = atoi(optarg);
        else if (c == 'R')
            initRad = atoi(optarg); // cooperativity


    }



    /////-------INIITIALIZE RNG------///

    const gsl_rng_type * T;
	gsl_rng * r;
	//---------Random Number Generator

	gsl_rng_env_setup();
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	//int sysRandom;
	gsl_rng_set(r, seedID); ///initialized rng
	//default_random_engine generator;
	//std::mt19937 gen( rd()) ;
	//auto gen = std::default_random_engine(seedID);

	/////-----INITIALIZE ADDITIONAL PARAMETERS_-----------
	vector <int> cellList;
	vector<std::tuple<int, int, int,  int ,int ,int> > lHistory;

	long double patchArr[nDemes][nDemes] = {{0}};
	long double popArr[nDemes][nDemes] = {{0}}; 


	//--------INITIALIZE DATA FILES -----
	ofstream flog, fprof,fpatch,fhist;
	time_t time_start;
	clock_t c_init = clock();
	struct tm * timeinfo;
	char buffer [80];


    time (&time_start);
	timeinfo = localtime (&time_start);


	strftime (buffer,80,"%F-%H-%M-%S",timeinfo);

	ostringstream date_time, Rstr, Istr;
	date_time << buffer;
	Rstr << initRad;
	Istr << seedID;
	string param_string =  "R"+Rstr.str()+"_I"+Istr.str()+"_";



	string logName = "log_" + param_string + date_time.str() + ".txt";
	string profName = "prof_" + param_string + date_time.str() + ".txt";
	string patchName = "patch_" + param_string + date_time.str() + ".txt";
	string historyName = "hist_" + param_string + date_time.str() + ".txt";



	//string folder = "sim_data/";
	string folder="";
    flog.open(folder+logName);
    fprof.open(folder+profName);
    fhist.open(folder+historyName);
    fpatch.open(folder+patchName);

    //int labelUnique =1;
    //int x = int(nDemes/2)+1;
    //int y =int(nDemes/2)+1;
	//popArr[x][y] =1;  
	//cellList.push_back(x*nDemes + y);

    int labelUnique =1;
	for(int i = 0; i < int(nDemes); i++){
		for(int j = 0; j < int(nDemes); j++){
			if ( round(sqrt( abs(i-int(nDemes*.5))*abs(i-int(nDemes*.5)) + abs(j-int(nDemes*.5))*abs(j-int(nDemes*.5))) ) < initRad )
			{
				popArr[i][j] =labelUnique ;  
				labelUnique+=1;


			}
		}
	}


	for(int i = 0; i < int(nDemes); i++){
		for(int j = 0; j < int(nDemes); j++){
			if ( round(sqrt( abs(i-int(nDemes*.5))*abs(i-int(nDemes*.5)) + abs(j-int(nDemes*.5))*abs(j-int(nDemes*.5))) ) < initRad )
			{
				if (checkEmpty(popArr, i,j,nDemes)==1){

					cellList.push_back(i*nDemes + j);

				}

			}
		}
	}


	/*for(int i = 0; i < int(nDemes); i++){
		for(int j = 0; j < int(nDemes); j++){

			patchArr[i][j] = 1 ;  



		}
	}


	for(int p=0; p<patchP;p++){

		int centerX = 	gsl_random_uniform_int(r,nDemes);
		int centerY = 	gsl_random_uniform_int(r,nDemes);


		for(int i = 0; i < int(nDemes); i++){
			for(int j = 0; j < int(nDemes); j++){
				if ( round(sqrt( abs(i-centerX )*abs(i-centerX ) + abs(j-centerY )*abs(j- centerY )) ) < patchS )
				{
					patchArr[i][j] = 1  ;  


				}
			}
		}
	}*/



	int linPatchD = int(nDemes/(2*patchS +patchP ));

	for(int xr = -int(linPatchD/2)+1; xr< int(linPatchD/2); xr++){
		for(int yr = -int(linPatchD/2)+1; yr< int(linPatchD/2); yr++){
			int xc= int(nDemes/2) +xr*(2*patchS +patchP);
			int yc= int(nDemes/2) +yr*(2*patchS +patchP);
			for(int i = 0; i < int(nDemes); i++){
				for(int j = 0; j < int(nDemes); j++){

					if ( round(sqrt( abs(i-int(xc))*abs(i-int(xc)) + abs(j-int(yc))*abs(j-int(yc))) ) < patchS)
					{
						patchArr[i][j]=1;
						


					}

				}


			}
		}
	}






	int i;
	int j;

	//cout << cellList.size()<<endl;
	for (int dt = 0 ; dt < nGen; dt++ ){
		cout << cellList.size()<<endl;


		//uniform_int_distribution<int> distribution_g(0, cellList.size() ) ;

		int cellID = cellList[ gsl_rng_uniform_int(r, cellList.size() ) ] ;



		j = cellID % nDemes;
		i =  int(cellID/nDemes);
		vector <int> neighborPick{0,0};
		neighborPick = pickEmpty(popArr, i,j,nDemes, r);





		popArr[neighborPick[0]][neighborPick[1]] = popArr[i][j];


		if ( patchArr[i][j] ==0 ){

			popArr[i][j] = 0;


			cellList.erase(remove(cellList.begin(), cellList.end(), i*nDemes+j), cellList.end()); 


			int arr[2] = {i,j}; // store current lattice coordinates
			int neighbDiff[4][2] = {{0,1},{0,-1},{1,0},{-1,0}}; //directions for negihbors
			int neighbCoor[4][2]; // initialize array for neighbor coordinates
			vector <int> neighbEmpty;
			//int fullCount =0;
			int ne =0;
			int neighbID; 
			for (int ne =0;ne<4;ne++){

				neighbCoor[ne][0] = (arr[0] + nDemes+neighbDiff[ne][0]) % nDemes; 
				neighbCoor[ne][1] = (arr[1] + nDemes+neighbDiff[ne][1]) % nDemes;
				neighbID = neighbCoor[ne][0] * nDemes + neighbCoor[ne][1];

				if (( popArr[neighbCoor[ne][0]][neighbCoor[ne][1]]  > 0) &&  (checkEmpty(popArr, neighbCoor[ne][0] ,neighbCoor[ne][1], nDemes) == 1 ) ){
					
					cellList.push_back(neighbID); 


				}

			}






		} else{
			lHistory.push_back(    make_tuple( popArr[i][j],i,j, neighborPick[0],neighborPick[1],  dt)   );
		}

		if ( checkEmpty(popArr, neighborPick[0],neighborPick[1], nDemes) == 1 && patchArr[i][j]==1 ) {


			cellList.push_back(neighborPick[0]*nDemes+ neighborPick[1] );


		}
	

		int arr1[2] = {neighborPick[0],neighborPick[1]}; // store current lattice coordinates
		int neighbDiff[4][2] = {{0,1},{0,-1},{1,0},{-1,0}}; //directions for negihbors
		int neighbCoor[4][2]; // initialize array for neighbor coordinates
		vector <int> neighbEmpty;
		//int fullCount =0;
		int ne =0;
		int neighbID; 

		for (int ne =0;ne<4;ne++){

			neighbCoor[ne][0] = (arr1[0] + nDemes+neighbDiff[ne][0]) % nDemes; 
			neighbCoor[ne][1] = (arr1[1] + nDemes+neighbDiff[ne][1]) % nDemes;
			neighbID = neighbCoor[ne][0] * nDemes + neighbCoor[ne][1];

			if (( find(cellList.begin(), cellList.end(), neighbID) != cellList.end() ) &&  (checkEmpty(popArr, neighbCoor[ne][0] ,neighbCoor[ne][1], nDemes) == 0 ) ){
				
				cellList.erase(remove(cellList.begin(), cellList.end(), neighbID), cellList.end()); 


			}

		}



	}	



	for(int i=0;i < nDemes; i++){
    	for(int j =0;j < nDemes; j++){


    		fprof << i << " " << j<< " " << popArr[i][j] << endl;
    		fpatch << i << " " << j<< " " << patchArr[i][j] << endl;
    		
		}

    }
    for (auto [ pID, pX, pY,cX,cY,T ] : lHistory)
    {
      fhist << pID << " " << pX << " " << pY << " "<<cX <<" "<<cY<< " "<< T << "\n"<<endl;
    }


	flog.close();
    fprof.close();


    clock_t c_fin = clock();
    double run_time = double(c_fin - c_init)/CLOCKS_PER_SEC;



    cout << "Finished!" << "\n";
    cout << seedID << endl;
    cout << "Finished in " << run_time << " seconds \n";

	puts (buffer);







}
