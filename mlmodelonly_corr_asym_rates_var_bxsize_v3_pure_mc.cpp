#include <iostream>
#include <fstream> // ofstream //
#include <iomanip> // setw(6) //
#include <cmath>
#include <string>
#include <stdlib.h>     // malloc, free, rand, atoi, atof //
#include <chrono>
#include <math.h>
#include <sstream>
//#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/random_device.hpp>
#include <boost/random/variate_generator.hpp>
#include <sys/time.h>
#include <vector>
//"/home/mayleonardmc/boost/boost_1_66_0/
//#include <boost/random/uniform_real.hpp>
// /home/shann87/boost/boost_1_66_0
// /home/shann87/boost/boost_1_66_0/stage/lib

//#include "/home/shann87/mayleonardmc/function/reaction_ml_only.cpp"   //call the function movementt()
//#include "/home/shann87/mayleonardmc/function/movement.cpp" //call the function reaction_lv()
//#include "/home/shann87/mayleonardmc/function/random_ini.cpp"  //call the function random_ini()
// #include<LV2D_onsites.h>

using namespace boost;
// using namespace std;


int time(){
struct timeval tp;
gettimeofday(&tp, NULL);
int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
return ms;
}

boost::random::mt19937 rng(time());
uniform_real<> dist(0.0,1.0);
variate_generator<mt19937,uniform_real<> > randomm(rng, dist);

inline double std_rand(){

    return randomm();
}

void random_ini(int* specc[], double density, int bx_y, int tot_bx_x ){   //generate a randomly distributed configuration of num paricle of type on a boxsize*boxsize 2D square lattice
	int position[2];
	int m=1;
	double randno;
	//std::cout << "Start Loop ini";
	if (density>0.33)
	 std::cout << "Density is wrongly entered";

    	for(int i=0; i<bx_y; i++ ){
		for(int j=0; j<tot_bx_x; j++){
			// sectioned initial conditions
			/*if (((i <= boxsize/3) && (j <= boxsize/3)) || ((i > boxsize/3 && i <= 2*boxsize/3) && (j > boxsize/3 && j <= 2*boxsize/3)) || ((i >= 2*boxsize/3) && (j >= 2*boxsize/3)))
			{
				specc[i][j] = 1;
			}
			else if	(((i <= boxsize/3) && (j > boxsize/3 && j <= 2*boxsize/3)) || ((i > boxsize/3 && i <= 2*boxsize/3) && (j >= 2*boxsize/3)) || ((i >= 2*boxsize/3) && (j <= boxsize/3)))
			{
				specc[i][j] = 2;
			}
			else if (((i <= boxsize/3) && (j >= 2*boxsize/3)) || ((i > boxsize/3 && i <= 2*boxsize/3) && (j <= boxsize/3)) || ((i >= 2*boxsize/3) && (j > boxsize/3 && j <= 2*boxsize/3)))
			{
				specc[i][j] = 3;
			}*/

				randno=(std_rand()+0.0)/m;
				if(randno >= 0 && randno < density)
					specc[i][j] = 1;
				else if (randno >= density && randno < 2*density)
					specc[i][j] = 2;
				else if (randno >= 2*density && randno < 3*density)
					specc[i][j] = 3;
				else
					specc[i][j] = 0;//}*/
		}
 	}

}

void stripe_ini(int* specc[], double density, int bx_y, int tot_bx_x ){   //generate a randomly distributed configuration of num paricle of type on a boxsize*boxsize 2D square lattice
	int position[2];
	int m=1;
	double randno;
	//std::cout << "Start Loop ini";
	if (density>0.33)
	 std::cout << "Density is wrongly entered";

    	for(int i=0; i<bx_y; i++ ){
		for(int j=0; j<tot_bx_x; j++){
			// sectioned initial conditions
			/*if (((i <= boxsize/3) && (j <= boxsize/3)) || ((i > boxsize/3 && i <= 2*boxsize/3) && (j > boxsize/3 && j <= 2*boxsize/3)) || ((i >= 2*boxsize/3) && (j >= 2*boxsize/3)))
			{
				specc[i][j] = 1;
			}
			else if	(((i <= boxsize/3) && (j > boxsize/3 && j <= 2*boxsize/3)) || ((i > boxsize/3 && i <= 2*boxsize/3) && (j >= 2*boxsize/3)) || ((i >= 2*boxsize/3) && (j <= boxsize/3)))
			{
				specc[i][j] = 2;
			}
			else if (((i <= boxsize/3) && (j >= 2*boxsize/3)) || ((i > boxsize/3 && i <= 2*boxsize/3) && (j <= boxsize/3)) || ((i >= 2*boxsize/3) && (j > boxsize/3 && j <= 2*boxsize/3)))
			{
				specc[i][j] = 3;
			}*/
			if (j < ceil(tot_bx_x/3))
				specc[i][j]=1;
			else if (j < ceil(2*tot_bx_x/3))
				specc[i][j]=2;
			else
				specc[i][j]=3;	
		}
 	}

}



/*void reaction_mod(int* specc[], int boxsize){   //generate a randomly distributed configuration of num paricle of type on a boxsize*boxsize 2D square lattice
	int position[2];
	int m=1;
	double randno;
    	for(int i=int(boxsize/4); i<(0.75)*boxsize; i++ ){
		for(int j=int(boxsize/4); j<(0.75)*boxsize; j++){
			//std::cout << "Array Loop";
			randno=(std_rand()+0.0)/m;
				if(randno <= 0.33 )
					specc[i][j] = 1;
				else if (randno <= 0.67)
					specc[i][j] = 1;
				else if (randno > 0.67)
					specc[i][j] = 1;

		}
 	}

}
*/



void movementt(int arr[], int n, int bx_y, int tot_bx_x)
{
	//std::cout << "positions" << arr[0] << "_" << arr[1] << std::endl;    
    
    float x = std_rand();
    if (x <= .25){ // right
        arr[1] = arr[1]+1;
        if(arr[1]>=tot_bx_x) arr[1]=0;
    }
    else if (x <= .50 && x > .25){ // down
        arr[0] = arr[0]+1;
        if(arr[0]>=bx_y) arr[0]=0;
    }
    else if (x <= .75 && x > .50){ // left
        arr[1] = arr[1]-1;
        if(arr[1]<0) arr[1]=tot_bx_x-1;
    }
    else{ //  up
        arr[0] = arr[0]-1;
        if(arr[0]<0) arr[0]=bx_y-1;
    }
}


void density(int *arr[], double &den_A, double &den_B, double &den_C,  int xbox, int ybox, int xbox0){

	den_A=0;
	den_B=0;
	den_C=0;
	// Comment for the density plots
	for(int j=0;j< ybox; j++){
		for(int i=xbox0; i < xbox0+xbox; i++){
			switch (arr[j][i]) {
				  case 1 :
					    // Process for test = 1
    					den_A++;
    				  break;

				  case 2 :
    					// Process for test = 5
    				        den_B++;
    				  break;

				  case 3 :
    					// Process for test = 5
    				        den_C++;
    				  break;

			}


		}
	}
	// Comment ends
	den_A=den_A/((xbox)*ybox);
	den_B=den_B/((xbox)*ybox);
	den_C=den_C/((xbox)*ybox);

}



void cor(int *arr[], double *cxy[], int colength, int xbox, int ybox, int xbox0){
	long int xy[10];
	int positionx;
	int start=xbox0;
	int end=xbox0+xbox;
	int counter=0;
	bool counter_flag=true;
	//std::cout << "Corr functions start !"; 
    //colength =xbox;
	for(int x=0; x < colength; x++){
		//xy=0;
		for (int q = 0; q < 10; q++) // Using for loop we are initializing
		{
    			xy[q] = 0;
		}
        counter=0;
        counter_flag=true;
		for(int i=0; i<ybox; i++){
            if(counter_flag==true) 
            	{counter=0;}
            for(int j=start; j<end; j++){

				positionx = j+x;
				if(positionx >= end)
                    continue;
                else {
                    if(counter_flag ==true) 
                    	counter++;
				//if(positionx >= boxsize) positionx =  positionx - boxsize;
                    switch(arr[i][j]){//+arr[positionx][j]){
                        case 1:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: xy[0]++; //aa
                                    break;
                                case 2: xy[1]++; //ab
                                    break;
                                case 3: xy[2]++; //ac
                                    break;
                                case 0:	xy[3]++; //a0
                                    break;}
                            break;
                        case 2:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: xy[1]++; //ab
                                    break;
                                case 2: xy[4]++; //bb
                                    break;
                                case 3: xy[5]++; //bc
                                    break;
                                case 0:	xy[6]++; //b0
                                    break;}
                            break;
                        case 3:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: xy[2]++; //ac
                                    break;
                                case 2: xy[5]++; //bb
                                    break;
                                case 3: xy[7]++; //cc
                                    break;
                                case 0:	xy[8]++; //c0
                                    break;}
                            break;
                        case 0:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: xy[3]++; //a0
                                    break;
                                case 2: xy[6]++; //b0
                                    break;
                                case 3: xy[8]++; //c0
                                    break;
                                case 0:	xy[9]++; //00
                                    break;}

						//break;
                    }
                }

			}
			// std::cout << "Counter: " << counter << std::endl;
			counter_flag=false;	
		}

		for (int l= 0; l< 10; l++){
			if (l ==0 || l ==4 || l ==7 || l ==9)
				//{
				{if (std::isnan(cxy[l][x]+(xy[l]+0.0)/((counter)*ybox)))
					std::cout << xy[l] << " l value " << l << " cxy " << cxy[l][x];
				cxy[l][x] = cxy[l][x]+(xy[l]+0.0)/((counter)*ybox);}
				 //std::cout << isnan(cxy[l][x]);}
			else
				cxy[l][x] = cxy[l][x]+(xy[l]+0.0)/((counter)*ybox);

		}
	}
//std::cout << "Corr functions work!"; 		
}
// ######################  CORR WITH AVG SUBSTRACTED ##########################################################################
void CumMinusAverage(int *arr[], double *cxy[], int colength, int xbox, int ybox, int xbox0){
	//  At the eleventh position we will have a[r], b[r] , c[r], \phi[r]
	// using namespace std;

	// std::cout << "Corr functions start !"<< std::flush; 

	long int xy[10+4];
	
	int positionx;
	int start=xbox0;
	int end=xbox0+xbox;
	int counter=0;
	bool counter_flag=true;
	//std::cout << "Corr functions start !"; 
    //colength =xbox;
	for(int x=0; x < colength; x++){
		//xy=0;
		for (int q = 0; q < (10+4); q++) // Using for loop we are initializing
		{
    			xy[q] = 0;
		}
        counter=0;
        counter_flag=true;
		for(int i=0; i<ybox; i++){
            if(counter_flag==true) 
            	{counter=0;}
            for(int j=start; j<end; j++){

				positionx = j+x;
				if(positionx >= end)
                    continue;
                else {
                    if(counter_flag ==true) 
                    	counter++;
				//if(positionx >= boxsize) positionx =  positionx - boxsize;
                    switch(arr[i][j]){//+arr[positionx][j]){
                        case 1:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: {xy[0]++; xy[10]++;//aa
                                    break;}
                                case 2: {xy[1]++; xy[11]++;//ab
                                    break;}
                                case 3: {xy[2]++; xy[12]++;//ac
                                    break;}
                                case 0:	{xy[3]++; xy[13]++;//a0
                                    break;}}
                            break;
                        case 2:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: {xy[1]++; xy[10]++;//ab
                                    break;}
                                case 2: {xy[4]++; xy[11]++;//bb
                                    break;}
                                case 3: {xy[5]++; xy[12]++;//bc
                                    break;}
                                case 0:	{xy[6]++; xy[13]++;//b0
                                    break;}}
                            break;
                        case 3:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: {xy[2]++; xy[10]++;//ac
                                    break;}
                                case 2: {xy[5]++; xy[11]++;//bb
                                    break;}
                                case 3: {xy[7]++; xy[12]++;//cc
                                    break;}
                                case 0:	{xy[8]++; xy[13]++;//c0
                                    break;}}
                            break;
                        case 0:
                            switch(arr[i][positionx]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                                case 1: {xy[3]++; xy[10]++;//a0
                                    break;}
                                case 2: {xy[6]++; xy[11]++;//b0
                                    break;}
                                case 3: {xy[8]++; xy[12]++;//c0
                                    break;}
                                case 0:	{xy[9]++; xy[13]++;//00
                                    break;}}

						//break;
                    }
                }

			}
			// std::cout << "Counter: " << counter << std::endl;
			counter_flag=false;	
		}
		
		for (int l= 0; l< (10+4); l++){
			if (l ==0 || l ==4 || l ==7 || l ==9)
				//{
				{if (std::isnan(cxy[l][x]+(xy[l]+0.0)/((counter)*ybox)))
					std::cout << xy[l] << " l value " << l << " cxy " << cxy[l][x];
				cxy[l][x] = cxy[l][x]+(xy[l]+0.0)/((counter)*ybox);}
				 //std::cout << isnan(cxy[l][x]);}
			else
				cxy[l][x] = cxy[l][x]+(xy[l]+0.0)/((counter)*ybox);

		}
	}
//std::cout << "Corr functions work!"; 		
}
// ######################  CORR WITH AVG SUBSTRACTED ##########################################################################
// ######################  CORR WITH AVG SUBSTRACTED ##########################################################################
void CumMinusAverageYdirection(int *arr[], double *cxy[], int colength, int xbox, int ybox, int xbox0){
	//  At the eleventh position we will have a[r], b[r] , c[r], \phi[r]
	// using namespace std;
	 
	long int xy[10+4];
	
	int positiony;
	int start=xbox0;
	int end=xbox0+xbox;
	int counter=0;
	bool counter_flag=true;
	
    //colength =xbox;

	for(int x=0; x < colength; x++){
		//xy=0;
		for (int q = 0; q < (10+4); q++) // Using for loop we are initializing
		{
    			xy[q] = 0;
		}
        counter=0;
        counter_flag=true;
		// for(int i=0; i<ybox; i++){
        for(int j=start; j<end; j++){
            if(counter_flag==true) 
            	{counter=0;}
            for(int i=0; i<ybox; i++){

				positiony = i+x;
				// cout << positiony << endl;
				if(positiony > (ybox-1)){
                    positiony=(ybox+positiony)%colength;
				}
				// if (positiony > (ybox-1)){
				// 	cout<<positiony<<endl;
				// 	cout << "x" << x << endl;
				// 	cout << "i" << i << endl;
				// }
				if(counter_flag ==true) 
                	counter++;
				//if(positionx >= boxsize) positionx =  positionx - boxsize;
                switch(arr[i][j]){//+arr[positionx][j]){
                    case 1:
                        switch(arr[positiony][j]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                            case 1: {xy[0]++; xy[10]++;//aa
                                break;}
                            case 2: {xy[1]++; xy[11]++;//ab
                                break;}
                            case 3: {xy[2]++; xy[12]++;//ac
                                break;}
                            case 0:	{xy[3]++; xy[13]++;//a0
                                break;}}
                        break;
                    case 2:
                        switch(arr[positiony][j]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                            case 1: {xy[1]++; xy[10]++;//ab
                                break;}
                            case 2: {xy[4]++; xy[11]++;//bb
                                break;}
                            case 3: {xy[5]++; xy[12]++;//bc
                                break;}
                            case 0:	{xy[6]++; xy[13]++;//b0
                                break;}}
                        break;
                    case 3:
                        switch(arr[positiony][j]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                            case 1: {xy[2]++; xy[10]++;//ac
                                break;}
                            case 2: {xy[5]++; xy[11]++;//bb
                                break;}
                            case 3: {xy[7]++; xy[12]++;//cc
                                break;}
                            case 0:	{xy[8]++; xy[13]++;//c0
                                break;}}
                        break;
                    case 0:
                        switch(arr[positiony][j]){//(speccc[i][j]!=0 && speccc[i][positionx]!=0) ab = ab + 1;
                            case 1: {xy[3]++; xy[10]++;//a0
                                break;}
                            case 2: {xy[6]++; xy[11]++;//b0
                                break;}
                            case 3: {xy[8]++; xy[12]++;//c0
                                break;}
                            case 0:	{xy[9]++; xy[13]++;//00
                                break;}}

					//break;
                }
            

			}
			// std::cout << "Counter: " << counter << std::endl;
			counter_flag=false;	
		}
		// cout << counter << "counter" << endl;
		for (int l= 0; l< (10+4); l++){
			// if (l ==0 || l ==4 || l ==7 || l ==9)
			// 	//{
			// 	{if (std::isnan(cxy[l][x]+(xy[l]+0.0)/((counter)*ybox)))
			// 		cout << xy[l] << " l value " << l << " cxy " << cxy[l][x] << "x" << x << endl;
			// 	cxy[l][x] = cxy[l][x]+(xy[l]+0.0)/((counter)*ybox);}
			// 	 //std::cout << isnan(cxy[l][x]);}
			// else
			// 	cxy[l][x] = cxy[l][x]+(xy[l]+0.0)/((counter)*ybox);
			// if (l ==0 || l ==4 || l ==7 || l ==9)
				if (std::isnan(cxy[l][x]+(xy[l]+0.0)/((counter)*ybox)))
					std::cout << xy[l] << " l value " << l << " cxy " << cxy[l][x] << "x" << x << std::endl;
				cxy[l][x] = cxy[l][x]+(xy[l]+0.0)/((counter)*ybox);
		}
	}
//std::cout << "Corr functions work!"; 		
}
// ######################  CORR WITH AVG SUBSTRACTED ##########################################################################

// 

void reaction_ml(int *specc[], double predation, double predRPS, double sex, double move, int* total_A, int* total_B, int* total_C, int* total_empty, int bx_y, int bx_x, int R, double asy_fact)
// void reaction_ml(int *specc[], double pred_A, double pred_B, double pred_C, double sex_A, double sex_B, double sex_C, double move_A, double move_B, double move_C, int* total_A, int* total_B, int* total_C, int* total_empty, int boxsize)				// this is there from when shannon wantted to look at diffrent predation and reproduction rates for A,B, and C
{
    //int where=0;
    int m=1;
    int boxsize = (1+R)*bx_x;
    int position[2];
    double total_rate;
    bool a_pred_rate= false;
    do{
        position[1] = (std_rand()+0.0)/m*boxsize;
        if(position[1] >= boxsize) position[1] = 0;
        position[0] = (std_rand()+0.0)/m*bx_y;
        if(position[0] >= bx_y) position[0] = 0;
    }while(specc[position[0]][position[1]] == 0);

	//std::cout << "positions" << position[0] << "_" << position[1] << std::endl;    
	//Half system simulate
	if (position[1] > (int)(bx_x-1))
		a_pred_rate = true;
	else
		a_pred_rate = false;
    // a_pred_rate=true;
	if (specc[position[0]][position[1]] == 1){
		// ASYMMETRIC RATES DEFINED HERE
		//std::cout << "where =1" << std::endl;
		// if (a_pred_rate==true){
		// 	 predation=0.1*predation;
		// }
        int position_new[2];
        position_new[0] = position[0];
        position_new[1] = position[1];
        movementt(position_new, 2, bx_y, boxsize);
		if ((specc[position_new[0]][position_new[1]] == 2) || (specc[position_new[0]][position_new[1]] == 3)){
			if (specc[position_new[0]][position_new[1]] == 2)
                if (a_pred_rate==true){
			        predation=asy_fact*predation;
		        }        
   //          total_rate=predation+move;
			// predation=predation/total_rate;
			//move=move/total_rate;
			}
			// else if (specc[position_new[0]][position_new[1]] == 0){
			// 	total_rate=sex+move;
			// 	sex=sex/total_rate;
			// 	}
		if (specc[position_new[0]][position_new[1]] == 2){
                    double rand_num = (std_rand()+0.0)/m;
		    if (rand_num <= predation/2){
				specc[position_new[0]][position_new[1]] = 0;
                        	(*total_B)--;
				(*total_empty)++;
			}
			else if ( rand_num > predation/2 && rand_num <= ((move+predation)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				}   		
		 //    else {
			//     int aaa = specc[position[0]][position[1]];
			//     specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		 //            specc[position_new[0]][position_new[1]] = aaa;
			// }
            // ASYMMETRIC RATES CHANGE BACK
		    predation=(1/asy_fact)*predation;
		}
        else if (specc[position_new[0]][position_new[1]] == 3){
		    double rand_num = (std_rand()+0.0)/m;
            if (rand_num <= predation/2){
			specc[position[0]][position[1]] = 0;
                        (*total_A)--;
			(*total_empty)++;
			}
		    else if ( rand_num > predation/2 && rand_num <= ((move+predation)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				} 
		 //    else {
			//     int aaa = specc[position[0]][position[1]];
			//     specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		 //            specc[position_new[0]][position_new[1]] = aaa;
			// }
		}
        else if (specc[position_new[0]][position_new[1]] == 0){
                    double rand_num = (std_rand()+0.0)/m;
            if (rand_num <= sex/2){
				specc[position_new[0]][position_new[1]] = 1;
                (*total_A)++;
				(*total_empty)--;
			}
			else if ( rand_num > sex/2 && rand_num <= ((sex+move)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				}   
		 //    else {
			//     int aaa = specc[position[0]][position[1]];
			//     specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		 //            specc[position_new[0]][position_new[1]] = aaa;
			// }
		}
		
	    }

    else if (specc[position[0]][position[1]] == 2){
    	//std::cout << "where =2" << std::endl;
		int position_new[2];
                position_new[0] = position[0];
                position_new[1] = position[1];
                movementt(position_new, 2, bx_y, boxsize);
		if ((specc[position_new[0]][position_new[1]] == 1) || (specc[position_new[0]][position_new[1]] == 3)){
			if (specc[position_new[0]][position_new[1]] == 1)
                if (a_pred_rate==true){
			        predation=asy_fact*predation;
		        }        
   //          total_rate=predation+move;
			// predation=predation/total_rate;
			//move=move/total_rate;
		}
		// else if (specc[position_new[0]][position_new[1]] == 0){
		// 	total_rate=sex+move;
		// 	sex=sex/total_rate;
		// 	}
        if (specc[position_new[0]][position_new[1]] == 3){
                    double rand_num = (std_rand()+0.0)/m;
		    if (rand_num <= predation/2){
			specc[position_new[0]][position_new[1]] = 0;
                        (*total_C)--;
			(*total_empty)++;
			}
		else if ( rand_num > predation/2 && rand_num <= ((move+predation)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				} 
		 // 	
		  //   else {
			 //    int aaa = specc[position[0]][position[1]];
			 //    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		  //           specc[position_new[0]][position_new[1]] = aaa;
				// }

		}
        else if (specc[position_new[0]][position_new[1]] == 1){
		    double rand_num = (std_rand()+0.0)/m;
            if (rand_num <= predation/2){
			specc[position[0]][position[1]] = 0;
                        (*total_B)--;
			(*total_empty)++;

			}
			else if ( rand_num > predation/2 && rand_num <= ((move+predation)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				} 
		 // 	


		  //   else {
			 //    int aaa = specc[position[0]][position[1]];
			 //    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		  //           specc[position_new[0]][position_new[1]] = aaa;
				// }
                // ASYMMETRIC RATES CHANGE BACK
		        predation=(1/asy_fact)*predation;
		}
        else if (specc[position_new[0]][position_new[1]] == 0){
                    double rand_num = (std_rand()+0.0)/m;
            if (rand_num <= sex/2){
			specc[position_new[0]][position_new[1]] = 2;
                        (*total_B)++;
			(*total_empty)--;
			}
			else if ( rand_num > sex/2 && rand_num <= ((move+sex)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				} 
		 // 	

		  //   else {
			 //    int aaa = specc[position[0]][position[1]];
			 //    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		  //           specc[position_new[0]][position_new[1]] = aaa;
				// }
		}
	    }


	else if (specc[position[0]][position[1]] == 3){
		//std::cout << "where =3" << std::endl;
                int position_new[2];
                position_new[0] = position[0];
                position_new[1] = position[1];
                movementt(position_new, 2, bx_y, boxsize);
		// if ((specc[position_new[0]][position_new[1]] == 1) || (specc[position_new[0]][position_new[1]] == 2)){
		// 	total_rate=predation+move;
		// 	predation=predation/total_rate;
		// 	//move=move/total_rate;
		// }
		// else if (specc[position_new[0]][position_new[1]] == 0){
		// 	total_rate=sex+move;
		// 	sex=sex/total_rate;
		// 	}
        if (specc[position_new[0]][position_new[1]] == 1){
                    double rand_num = (std_rand()+0.0)/m;
		    if (rand_num <= predation/2){
			specc[position_new[0]][position_new[1]] = 0;
                        (*total_A)--;
			(*total_empty)++;
			}
			else if ( rand_num > predation/2 && rand_num <= ((move+predation)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				} 
		 // 
		  //   else {
			 //    int aaa = specc[position[0]][position[1]];
			 //    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		  //           specc[position_new[0]][position_new[1]] = aaa;
				// }

		}
        else if (specc[position_new[0]][position_new[1]] == 2){
		    double rand_num = (std_rand()+0.0)/m;
            if (rand_num <= predation/2){
			specc[position[0]][position[1]] = 0;
                        (*total_C)--;
			(*total_empty)++;
			}
			else if ( rand_num > predation/2 && rand_num <= ((move+predation)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				} 
		 // 

		  //   else {
			 //    int aaa = specc[position[0]][position[1]];
			 //    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		  //           specc[position_new[0]][position_new[1]] = aaa;
				// }
		}
        else if (specc[position_new[0]][position_new[1]] == 0){
            double rand_num = (std_rand()+0.0)/m;
            if (rand_num <= sex/2){
				specc[position_new[0]][position_new[1]] = 3;
                (*total_C)++;
				(*total_empty)--;
			}
		    else if ( rand_num > sex/2 && rand_num <= ((move+sex)/2)){
			    int aaa = specc[position[0]][position[1]];
			    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		            specc[position_new[0]][position_new[1]] = aaa;				
				} 
		 // 

		  //   else {
			 //    int aaa = specc[position[0]][position[1]];
			 //    specc[position[0]][position[1]] = specc[position_new[0]][position_new[1]];
		  //           specc[position_new[0]][position_new[1]] = aaa;
				// }
		}

	    }
}

int main(int argc, char *argv[]){
    using namespace std;

    // Elementary definitions
    int bx_y, seed, corrsize;
    double density_A, density_B, density_C;
    seed = atoi(argv[1]);
    bx_y = atoi(argv[2]);
    density_A = density_B = density_C = atof(argv[3]);
    double predation, sex, move, predRPS, predvar, oldpred;
    predation = atof(argv[4]);
    sex = atof(argv[5]);
    move = atof(argv[6]);
    predRPS=atof(argv[7]);
    long int printtime =atoi(argv[8]);
    long int endtime =atoi(argv[9]);

       // CONTROL VARIABLES
    bool pred_mod_flag = true;
    long int control_time,quench_counter;
    //double osc_freq=141.6;
    int time_spawn_site=atoi(argv[13]); // This is also the control time variable
    int site_size = atoi(argv[14]);  // Size of the  driving sublattice
    int species_spawned_type=0;
    int number_site_loc=atoi(argv[15]);

    int p_latt_time=atoi(argv[11]);
    int p_corr_time=atoi(argv[12]);

    double R=atoi(argv[17]);
    int bx_x =atoi(argv[16]);
    double asy_fact=atof(argv[18]);
    //cout << "variable read correct" << endl; 
    long int k=0;
    int counter =0;
    //double invboxsizesq =(1.0/(bx_y*bx_x)); // singularity possible

    long long int x = 0;
    long int intermediate_ppm_steps =0;
    int total_A, total_B, total_C, total_num, total_empty, tot_stab_population;
    total_A = bx_y * (1+R)*bx_x * density_A;
    total_B = total_C = total_A;
    tot_stab_population=0;

    total_num = total_A + total_B + total_C;
    total_empty = (bx_y * (1+R)*bx_x) - total_num;
    


    //cout << "Test 01\n";
    double t = 0.000000;
    long int timestep = 0;
    int position[2];
    int position_new[2];
    int timer=0;

    //



    // compute how many steps to look at from endtime and print time
    //intermediate_ppm_steps =pow(10,(floor(log10(endtime-printtime)-1)));
    //srand(seed);
    //cout << intermediate_ppm_steps << endl;
    int *specc[bx_y];
    for(int i=0; i<bx_y; i++){
        specc[i] = new int[(int)(1+R)*bx_x];
    }
    for(int i=0; i<bx_y; i++){
        for(int j=0; j<(1+R)*bx_x; j++){
            specc[i][j] = 0;
        }}
    // corr array
    /*    
    int *cspecc[bx_y];
    for(int i=0; i<bx_y; i++){
      		cspecc[i] = new int[bx_x];
    }
    for(int i=0; i<bx_y; i++){
        for(int j=0; j<bx_x; j++){
            cspecc[i][j] = 0;
        }}
    */    
    
    //double *time_corr;
    //for(int i=0; i<(endtime-printtime);i++){
	//time_corr = new double[endtime-printtime];


	double den_A, den_B, den_C;
	den_A = den_B = den_C =0;
	cout << den_A << den_B << den_C;

	int corrlength_sym=(int)(bx_x);
	int corrlength_asym=(int)(bx_x*R);
	
	double *scxy[14];
	for(int i=0; i<14; i++){
       		scxy[i] = new double[corrlength_sym];
    	}

    double *acxy[14];
	for(int i=0; i<14; i++){
       		acxy[i] = new double[corrlength_asym];
    	}
	

	int *new_specc[bx_y];
		for(int i=0; i<bx_y; i++){
        		new_specc[i] = new int[(int)(1+R)*bx_x];
    		}

	
    // Flags --
    bool lattice_print_flag = true;
    bool initial_time_flag = false;
    bool lattice_mod_flag= false;
    bool quench_flag=false;
    bool space_corr_flag =true;
    //bool quench_con=true;
    //


	// Extinction flags
	bool extinction_reached=false;

// Define using random numbers the location of the sites where we have asymmetrical rates
// Test one system where we just divide the system into two and implement asymmetric rates
// in one half.


/*

    // Implement 2d system in certain regions randomly.
    // Here the site size gives a square on which these reactions are to be performed
    //
    //int site_size=8;
    bool implement_a_rates=true;
    //reaction_lv_2d();
    bool *asy_latt[bx_y];
    if(implement_a_rates == true){
       // bool *asy_latt[bx_y];
        for(int i=0; i<bx_y; i++){
            asy_latt[i] = new bool[(1+R)*bx_x];
        }
        for(int i=0; i<bx_y; i++){
            for(int j=0; j<(1+R)*bx_x; j++){
                asy_latt[i][j] = false;
        }}
        int position[2];

        for(int l=0; l<number_site_loc;l++){
           position[0]=(std_rand()+0.0)*bx_y;
           position[1]=(std_rand()+0.0)*(1+R)*bx_x;
                for(int i=position[0];i<position[0]+site_size;i++)
                    for(int j=position[0];j<position[0]+site_size;j++)
                        asy_latt[i][j]==true;
        }

    }
	

*/

    //cout << " Array dim read correct" << endl;
    long int species_count =0;
    random_ini(specc, density_A, bx_y, (1+R)*bx_x);
	// stripe_ini(specc, density_A, bx_y, (1+R)*bx_x);

//   correlation file
/*
     stringstream ss2;
     ss2 << argv[10] << "S" << "_" << bx_y << "X" << (1+R)*bx_x << "_" << predation << "_"  << sex << "_" << move << "corrtime" <<  ".txt";
     string name2;
     name2 = ss2.str();
     fstream data2(name2.c_str(), ofstream::out | ofstream::app | ofstream::in);
*/
//

//	double random;
//	oldpred=predation;
//	predvar=0.1*predation;





   while(timestep <= endtime){
        total_num = total_A + total_B + total_C;
        t = t + 1.00000/total_num;
		species_count=0;


		for(int q=0; q<(bx_y*(1+R)*bx_x); q++){
        		reaction_ml(specc, predation, predRPS, sex, move, &total_A, &total_B, &total_C, &total_empty, bx_y, bx_x, R, asy_fact);
        		//cout << "reaction ml done correct" << q << "_" << timestep << endl;
			/*if (timestep % time_spawn_site == 0){
				for (int i=(int)(boxsize/2-site_size/2); i < (int)(boxsize/2+site_size/2) +1 ; i++){
					for (int j=(int)(boxsize/2-site_size/2); j < (int)(boxsize/2+site_size/2) +1 ; j++)
						specc[i][j] = species_spawned_type+1;
				}
				species_spawned_type =(species_spawned_type +1 ) % 3 ;
			}*/
    		}



		//if(initial_time_flag ==true)
		//	cout<< initial_time_flag;

		/*if((timestep > printtime) && initial_time_flag ==false){
			//initial_time_flag =true;
			tot_stab_population = total_A + total_B + total_C;
			control_time=int(osc_freq*tot_stab_population);
			cout << timestep << "control time " << control_time << endl;
			cout << "tot pop" << tot_stab_population << endl;
			//quench_flag = true;
			//cout << "Condition satis" << endl;
		}*/
		//timer=time();
		//cout<< timer << endl;
		/*if((timestep % control_time == 0) && (timestep > printtime)){
				//cout <<  timestep << "timestep";
				random=0;
				counter++;
				lattice_mod_flag =true;
				//reaction_mod(specc, boxsize);
				//cout << "Counter" << counter << endl;
				for(int u=64; u<192; u++){
					//cout <<  timestep << "timestep";
        				for(int v=64; v<192; v++){
						//specc[u][v] = 2;
						//cout <<  timestep << "timestep";
						random=(std_rand()+0.0);
						//random++;
						//if(random < 0.33)
						  if(random < (total_A*invboxsizesq))
							specc[u][v] = 1;
						else if((random >= (total_A*invboxsizesq)) && (random < ((total_A+total_B)*invboxsizesq)))
							specc[u][v] = 2;
						else if((random >= ((total_A+total_B)*invboxsizesq)) && (random < ((total_A+total_B+total_C)*invboxsizesq)))
							specc[u][v] = 3;
						else
							specc[u][v] = 0;

				 }}


		}*/
		/*
		quench_counter++;

		if((timestep % cosc_freq == 0)){
			if(quench_flag == true){
				quench_counter=0;
				quench_flag=false;
				predation=predvar;
			}
			else if(quench_flag==false){
				quench_counter=0;
				quench_flag=true;
				predation=oldpred;
			}

		}
		*/


	// SPACE CORR CODE

	if (space_corr_flag == true && timestep % p_corr_time ==0){

		for(int i=0; i<14; i++){
			// This array should have the combined correlation size of the symmetric and asymmetric region
        		for(int j=0; j<(corrlength_sym); j++){
            			scxy[i][j] = 0;
        	}}
//void cor(int *arr[], double *cxy[], int colength, int xbox, int ybox, int xbox0){
		// cor(specc, scxy, corrlength_sym, bx_x, bx_y, 0);

		//	changes for a(r)  and  y direction
		// cout << "entering cum function" << endl;
		CumMinusAverage(specc, scxy, corrlength_sym, bx_x, bx_y, 0);
		// cout << "exiting cum funciton" << endl;
		

		CumMinusAverageYdirection(specc, scxy, (int)(bx_y/2), bx_x, bx_y, 0);
		//  Divide correlation function bu two
		for(int i=0;i<(int)(bx_y/2);i++)
		{
			for(int j=0;j<14;j++)
				scxy[j][i]=scxy[j][i]/2;
		}

		//  changes for a(r)  and  y direction

		//cout << "Corr arrays correct" << endl;
		for(int i=0; i<14; i++){
			// This array should have the combined correlation size of the symmetric and asymmetric region
        		for(int j=0; j<(corrlength_asym); j++){
            			acxy[i][j] = 0;
        	}}

		// cor(specc, acxy, corrlength_asym,R*bx_x, bx_y, bx_x);
	
		//	changes for a(r)  and  y direction
		CumMinusAverage(specc, acxy, corrlength_asym,R*bx_x, bx_y, bx_x);
		CumMinusAverageYdirection(specc, acxy, (int)(bx_y/2),R*bx_x, bx_y, bx_x);
		//  Divide correlation function bu two
		for(int i=0;i<(int)(bx_y/2);i++)
		{
			for(int j=0;j<14;j++)
				acxy[j][i]=acxy[j][i]/2;
		}
		// //  changes for a(r)  and  y direction
		

		/*for(int i=0;i < boxsize; i++){
			for(int j=0;j < boxsize; j++){
				new_specc[i][j] = specc[j][i];
			}
		}*/
		
        double den_A1,den_B1,den_C1;
        double den_A2,den_B2,den_C2;
		density(specc, den_A1, den_B1, den_C1 ,bx_x, bx_y, 0);

		density(specc, den_A2, den_B2, den_C2, R*bx_x, bx_y, bx_x);

		if (den_A1 ==0 || den_B1 ==0 || den_C1 ==0 || den_A2 ==0 || den_B2 ==0 ||den_C2 ==0 )
			extinction_reached=true;


		//cout << "den vars work correct" << endl;

		stringstream ssC;
        		ssC << argv[10] << "S_CORR_DEN" << "_" << bx_y << "X" << (1+R)*bx_x << "_" << predation << "_"  << sex << "_" << move << "_RPS"<< predRPS << "_" << setfill('0') << setw(9) << timestep << ".txt";
        		string nameC;
        		nameC = ssC.str();
        		fstream dataC(nameC.c_str(), ofstream::out | ofstream::app | ofstream::in);

			dataC << den_A1 << "," <<  den_B1 << "," << den_C1 << "," << den_A2 << "," <<  den_B2 << "," << den_C2 << "," << endl;

        for(int i=0;i < corrlength_sym; i++){
			for(int l=0; l < 14 ; l++){
				dataC << (scxy[l][i]) << " ";
			}
			dataC << endl;
		}
		for(int i=0;i < corrlength_asym; i++){
			for(int l=0; l < 14 ; l++){
				dataC << (acxy[l][i]) << " ";
			}
			dataC << endl;
		}
		dataC.close();
		//out << "corr_den printed correct" << endl;	

		//	changes for a(r)  and  y direction
		stringstream ssX;
		ssX << argv[10] << "S_CumMinAvg" << "_" << bx_y << "X" << (1+R)*bx_x << "_" << predation << "_"  << sex << "_" << move << "_RPS"<< predRPS << "_" << setfill('0') << setw(9) << timestep << ".txt";
		string nameX;
		nameX = ssX.str();
		fstream dataX(nameX.c_str(), ofstream::out | ofstream::app | ofstream::in);

		// dataX << den_A1 << "," <<  den_B1 << "," << den_C1 << "," << den_A2 << "," <<  den_B2 << "," << den_C2 << "," << endl;

        for(int i=0;i < corrlength_sym; i++){
			for(int l=10; l < 14 ; l++){
				dataX << (scxy[l][i]) << " ";
			}
			dataX << endl;
		}
		for(int i=0;i < corrlength_asym; i++){
			for(int l=10; l < 14 ; l++){
				dataX << (acxy[l][i]) << " ";
			}
			dataX << endl;
		}
		dataX.close();

	}

// 		//	changes for a(r)  and  y direction	

// 	}


	// SPACE CORR CODE


	// correlation changes
/*	if(timestep > printtime){

		//long int counter =0;
		//double time_count =0;

		if(initial_time_flag ==false){
			tot_stab_population = total_A + total_B + total_C;
			corrsize=int((endtime));
			//corrsize=int((endtime)/tot_stab_population);
			//corrsize=int((endtime-printtime)/tot_stab_population)+10;
			time_corr = new double[corrsize];
			initial_time_flag =true;
			for(int i=0; i<boxsize; i++){
        			for(int j=0; j<boxsize; j++){
            				cspecc[i][j] = specc[i][j];
        		}}
			time_corr[0]=1.0;
			cout << time_corr[0]<< "This should be 1";
		}
		//else if(timestep % (tot_stab_population) == 0){
		else {
			for(int i=0; i<boxsize; i++){
        			for(int j=0; j<boxsize; j++){
            				if(cspecc[i][j] == specc[i][j]){
						time_corr[k]++;
					}
	        		}
			}
			time_corr[k]= time_corr[k]/(float(boxsize*boxsize));
		}

			//if(timestep % (tot_stab_population) == 0){
                 	data2 << k << " " << time_corr[k] << " " << (total_A*invboxsizesq) << " " << (total_B*invboxsizesq) << " " << (total_C*invboxsizesq) << endl;
			k++;
				//}




	}
*/






    //     // correlation changes end
    //    if(timestep >= 000 && lattice_print_flag == true){
	//    if((timestep % p_latt_time ==0)||lattice_mod_flag==true){	//comment this block of code in order to measure the population density of species A
    //                     //           /home/shann87/mayleonardmc/testdata/

	// 		lattice_mod_flag=false;
    //         		stringstream ss;
    //         		ss << argv[10] << "S" << "_" << bx_y << "X" << (1+R)*bx_x << "_" << predation << "_"  << sex << "_" << move << "_RPS"<< predRPS << "_" << setfill('0') << setw(9) << timestep << ".ppm";
    //         		string name;
    //         		name = ss.str();
    //         		fstream data(name.c_str(), ofstream::out | ofstream::app | ofstream::in);

    //         		data << "P3    " <<  (int)((1+R)*bx_x) << "   " << bx_y  << " 1" << endl;
    //         		for(int j=0; j<bx_y; j++){
    //            			for(int k=0; k<(1+R)*bx_x; k++){
    //              	  	if(specc[j][k] == 0){
    //                			data << 1 << " " << 1 << " " << 1 << "  ";
    //               	 	}
    //               	 	else if(specc[j][k] == 1){
    //                			data << 1 << " " << 0 << " " << 0 <<  "  ";
    //               	 	}
    //                		else if(specc[j][k] == 2){
    //                  			data << 0 << " " << 1 << " " << 0 << "  ";
    //               		 }
    //               		else if(specc[j][k] == 3){
    //                	  		data << 0 << " " << 0 << " " << 1 << "  ";
    //              		  }
    //           	  	 	}
    //           	 	 	data << endl;

    //        	  	}

	// 		data.close();


    //        }
	// }

	timestep++;


}
//}

        /*
			stringstream ss1;
            		ss1 << argv[10] << "S" << "_" << boxsize << "_" << predation << "_"  << sex << "_" << move << "_RPS"<< predRPS << ".txt";
            		string name1;
            		name1 = ss1.str();
            		fstream data1(name1.c_str(), ofstream::out | ofstream::app | ofstream::in);


            		for(int j=0; j<boxsize; j++){
               			for(int k=0; k<boxsize; k++){
                 	  	 	data1 << specc[j][k] << " ";
              	  	 	}
              	 	 data1 << endl;
			}


    			data1.close();
			/*
			stringstream ss2;
            		ss2 << argv[10] << "S" << "_" << boxsize << "_" << predation << "_"  << sex << "_" << move << "corrtime" <<  ".txt";
            		string name2;
            		name2 = ss2.str();
            		fstream data2(name2.c_str(), ofstream::out | ofstream::app | ofstream::in);

            		//data << "P3    " <<  boxsize << "   " <<  boxsize << " 1" << endl;

               		for(int k=0; k < (endtime-printtime); k++){
				if(timestep % 10 == 0){
                 	  	 	data2 << timestep << time_corr[k] << endl;
				}
			}
			*/
              	// 	data2.close();

            // stringstream ssdat;
            // ssdat << argv[10] << "Sdatafile" <<  ".txt";
            // //  ssdat << "datafile.txt";
            // string namedat;
            // namedat = ssdat.str();
            // fstream dataD(namedat.c_str(), ofstream::out | ofstream::app | ofstream::in);
            // dataD << "directory," << argv[10] << endl;
            // dataD << "boxsize in ydir," << bx_y << endl;
            // dataD << "boxsize in xdir i.e (1+R)*bx_x," << (1+R)*bx_x << endl;
            // dataD << "predation," << predation << endl;
            // dataD << "sex," << sex << endl;
            // dataD << "move," << move << endl;
            // dataD << "printtime," << printtime << endl;
            // dataD << "endtime," << endtime << endl;
            // dataD << "spawn_every," << time_spawn_site << endl;
            // dataD << "site_size," << site_size << endl;
            // dataD << "number_site_loc," << number_site_loc << endl;
            // dataD << "p_latt_time," << p_latt_time << endl;
            // dataD << "p_corr_time," << p_corr_time << endl;
            // dataD << "R, " << R << endl;

            // dataD.close();




	stringstream ssdat;
	ssdat << argv[10] << "Sdatafile" <<  ".txt";
	//  ssdat << "datafile.txt";
	string namedat;
	namedat = ssdat.str();
	fstream dataD(namedat.c_str(), ofstream::out | ofstream::app | ofstream::in);
	dataD << extinction_reached << endl;

	dataD << "directory," << argv[10] << endl;
	dataD << "boxsize in ydir," << bx_y << endl;
	dataD << "boxsize in xdir i.e (1+R)*bx_x," << (1+R)*bx_x << endl;
	dataD << "predation," << predation << endl;
	dataD << "sex," << sex << endl;
	dataD << "move," << move << endl;
	dataD << "printtime," << printtime << endl;
	dataD << "endtime," << endtime << endl;
	dataD << "spawn_every," << time_spawn_site << endl;
	dataD << "site_size," << site_size << endl;
	dataD << "number_site_loc," << number_site_loc << endl;
	dataD << "p_latt_time," << p_latt_time << endl;
	dataD << "p_corr_time," << p_corr_time << endl;
	dataD << "R, " << R << endl;
	dataD << "asy_fact, " << asy_fact << endl;

	dataD.close();

}
