#include "ibfs.h"
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define USE_GPU_ID		4 
int main(int args, char *argv[]) {
	
	MPI_Init(&args,&argv);
	int world_sz;
	MPI_Comm_size(MPI_COMM_WORLD, &world_sz);
	
	int my_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	const index_t gpu_id 	= USE_GPU_ID;
	ibfs(world_sz, my_id, gpu_id, args, argv);

	MPI_Finalize();
	return 0;
}

