#include <mpi.h>
#include "ibfs.h"
#include <iostream>

int reporter(double tm, int my_id, index_t iter)
{
	MPI_Barrier(MPI_COMM_WORLD);
	double max_tm;
	MPI_Reduce(&tm, &max_tm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(my_id==0)std::cout<<"Traversal-iter-"<<iter<<": "<<max_tm<<" second(s)\n";
	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}
