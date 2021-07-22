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
const int nDemes = 200; /// side length of square simulation lattice
int initRad = 10; // initial innoculation radius in demes
int patchS = 10;
int patchP = 30;
unsigned int seedID = 2; //rng sseed
unsigned int nGen = 25;




long double popArr[nDemes][nDemes] = {{0}}; ///array for pop count for each species at ecah lattice site
long double patchArr[nDemes][nDemes] = {{0}}; //holder array for pop count for each species at ecah lattice site before migration

//int patchCoor[patchP][2] = {{0}};




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
	default_random_engine generator;
	//std::mt19937 gen( rd()) ;
	auto gen = std::default_random_engine(seedID);

	/////-----INITIALIZE ADDITIONAL PARAMETERS_-----------
	vector <int> iRand;
	vector <int> jRand;
	vector<std::tuple<int, int, int,  int ,int ,int> > lHistory;  //parent ID, parent x, perent y, child x,child y,t




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



    //------INITIALIZE RADIAL INNOCULATION---------
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

	//for (int p=0; p<patchP; p++ ){
	//	patchCoor[p][0] = 
	//	patch



	///}
	////----initiaize PATTCHES ----///
	uniform_int_distribution<int> distribution_p(0,nDemes);
	for(int p=0; p<patchP;p++){

		int centerX = distribution_p(generator);
		int centerY = distribution_p(generator);


		for(int i = 0; i < int(nDemes); i++){
			for(int j = 0; j < int(nDemes); j++){
				if ( round(sqrt( abs(i-centerX )*abs(i-centerX ) + abs(j-centerY )*abs(j- centerY )) ) < patchS )
				{
					patchArr[i][j] = 1  ;  
					labelUnique+=1;


				}
			}
		}
	}







	//initialize random array for evoltuoion
	for(int x=0; x< nDemes; x++){

		jRand.push_back(x);
		iRand.push_back(x);
	}


	/////------MAIN LOOP-----
	int ii;
	int jj;
	

	for (int dt = 0 ; dt < nGen; dt++ ){

		shuffle(iRand.begin(),iRand.end(),gen);
		shuffle(jRand.begin(),jRand.end(),gen );

		for(int i =0;i<nDemes;i++){

			ii = iRand[i];
			for (int j=0;j<nDemes;j++){

				jj=jRand[j];

				if (popArr[ii][jj] >0){

					int arr[2] = {ii, jj}; // store current lattice coordinates
					int neighbDiff[4][2] = {{0,1},{0,-1},{1,0},{-1,0}}; //directions for negihbors
					int neighbCorr[4][2]; // initialize array for neighbor coordinates
					vector <int> neighbEmpty;
					

					int emptCount=0;
					for(int ne=0; ne <4; ne++){
						 ///find x,u coordinate of neighbor respecting periodic boundaries
						neighbCorr[ne][0] = (arr[0] + nDemes+neighbDiff[ne][0]) % nDemes; 
						neighbCorr[ne][1] = (arr[1] + nDemes+neighbDiff[ne][1]) % nDemes;

						if (popArr[ neighbCorr[ne][0] ][ neighbCorr[ne][1] ] ==0){

							neighbEmpty.push_back(ne);
							emptCount+=1;
						}

					}

					if (emptCount>0){
					
						uniform_int_distribution<int> distribution_e(0,emptCount-1);
						int nePick = distribution_e(generator);

						popArr[ neighbCorr[ neighbEmpty[nePick] ] [0] ] [ neighbCorr[ neighbEmpty[nePick] ] [1] ] =popArr[ii][jj];

						if (patchArr[ii][jj] ==0){

							popArr[ii][jj] = 0;
						} else{
							lHistory.push_back(    make_tuple( popArr[ii][jj],ii,jj, neighbCorr[ neighbEmpty[nePick] ] [0],neighbCorr[ neighbEmpty[nePick] ] [1],  dt)   );
						}
					}
				}
			
			}

		}
		
	}




	//for(int i=0; i < lHistory.size();i++){
	//	for( int j =0 j <)

	//		fhistory << lHistory[i] << endl;
	//}


    for (auto [ pID, pX, pY,cX,cY,T ] : lHistory)
    {
      fhist << pID << " " << pX << " " << pY << " "<<cX <<" "<<cY<< " "<< T << "\n"<<endl;
    }


	for(int i=0;i < nDemes; i++){
    	for(int j =0;j < nDemes; j++){

    		fprof << i << " " << j<< " " << popArr[i][j] << endl;
    		fpatch<< i << " " << j << " " << patchArr[i][j] << endl;
 
		}

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
