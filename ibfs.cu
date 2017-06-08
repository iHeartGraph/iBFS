//Graph format: Json based format
//Storage format: 
//struct{
//		int: src_ver
//		Arr: [ver_0|ver_1|ver_2|...]
//		Int: num_conn_ver
//	}
#include "graph.h"
#include "graph_reader.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "ibfs.h"

//template<typename index_t>
//void proc_desc(
//	std::string desc_str,
//	std::string &kron_addr,
//	index_t	&sw_level,
//	index_t &vert_count,
//	index_t	&edge_count
//){
//	//proc desc file
//	std::ifstream descfile(desc_str.c_str());
//	std::stringstream ss;
//	std::string line;
//	index_t	count=0;
//
//	if(descfile.is_open()){
//		while(std::getline(descfile, line)){
//			count++;
//			switch(count){
//				case 5:
//					ss.str("");
//					ss.clear();
//					ss<<line;
//					ss>>vert_count;ss>>edge_count;break;
//				case 11:
//					ss.str("");
//					ss.clear();
//					ss<<line;
//					ss>>kron_addr;ss>>sw_level;break;
//				default:
//					break;
//			}
//		}
//	}else{std::cout<<"File desc.dat missing\n"; exit(-1);}
//	std::cout<<"gpu vert count: "<<vert_count<<"\t"
//					 <<"gpu edge count: "<<edge_count<<"\n";
//}

int ibfs(int world_sz, int my_id, int gpu_id, int args, char *argv[]) {
	//typedef int vertex_t;
	//typedef int index_t;
	//typedef 	unsigned char 	depth_t;
	//typedef 	uint4 comp_t;
	
	std::cout<<"Input: /path/to/exe /path/to/beg "
				<<"/path/to/csr concurrent-count " 
				<<"switch-level bfs_count\n";
	if(args != 6){std::cout<<"Wrong input\n";exit(-1);}
	
	const char *beg_file=argv[1];
	const char *csr_file=argv[2];
	const int concurrent_count=atoi(argv[3]);
	const depth_t sw_level=atoi(argv[4]);
	const index_t src_count=atoi(argv[5]);
	graph_reader *g=new graph_reader(beg_file, csr_file, src_count);

	cudaSetDeviceFlags(cudaDeviceMapHost);
	const index_t sml_shed 	= 32; 
	const index_t lrg_shed 	= 512;
	index_t gpu_count		= 1;

	int num_groups = src_count / concurrent_count;
	int num_agg_bfs = concurrent_count;
	
	graph<vertex_t, index_t, depth_t> *graph_d 
		= new graph<vertex_t, index_t, depth_t>
				(
				 	g,
					sw_level,
					g->vert_count,
					g->edge_count,
					num_groups,
					num_agg_bfs,
					world_sz, my_id,
					gpu_id,
					sml_shed, 
					lrg_shed);

//	index_t step_sz = graph_d->vert_count/(THDS_NUM*BLKS_NUM);
//	if((step_sz%16 != 0)||(graph_d->vert_count%(THDS_NUM*BLKS_NUM))){
//		std::cout<<"Graph vertices sz fails " 
//						 <<"to be exactly divided by thread count\n";
//		return -1;
//	}
	
	std::cout<<"Alloc GPU space\n";
	graph_d->alloc_array();
	std::cout<<"In gpu bfs\n";
	graph_d->bfs_gpu_coalescing_mem();
//	graph_d->write_result();
	return 0;
}
