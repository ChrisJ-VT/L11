#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

int main(int argc, char **argv) {

	//MPI Initialization
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


	//need runnig tallies
	long long int Ntotal;
	long long int Ncircle;

	//seed random number generator
	double seed = rank;
	srand48(seed);

	for (long long int n = 0; n<1000000;n++){
		//generates tow random numbers
		double rand1 = drand48(); //returns a number between 1 and 0
		double rand2 = drand48();

		double x = -1 + 2*rand1; //shift to  [-1,1]
		double y = -1 + 2*rand2;
	
		if (sqrt(x*x+y*y) <= 1) Ncircle++;
		Ntotal++;	
	}

	double pi = 4.0*Ncircle/ (double) Ntotal;
	printf("Our estimate of pi is %f \n", pi);

	MPI_Finalize();

	return 0;
}