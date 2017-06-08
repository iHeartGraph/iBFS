#include "comm.h"

template <typename vertex_t, typename index_t, typename depth_t>
int graph<vertex_t, index_t, depth_t>::alloc_array()
{
	cudaSetDevice(gpu_id);	
	//long cpu_bytes= 0;
	long gpu_bytes= 0;
	
	const size_t cat_index_sz 	= sizeof(index_t)*BLKS_NUM*THDS_NUM;
	//const size_t blk_sz 				= sizeof(index_t)*BLKS_NUM*num_agg_bfs;
	const size_t merge_depth_sz	= sizeof(depth_t)*vert_count*num_agg_bfs;
	const size_t edge_sz 				= sizeof(vertex_t)*edge_count;
	const size_t vert_sz 				= sizeof(vertex_t)*vert_count;
	const size_t index_sz 			= sizeof(index_t)*vert_count;
	const size_t comp_sz 				= sizeof(comp_t)*vert_count;
	const size_t src_sz					= sizeof(vertex_t)*src_count;
	gpu_bytes+= edge_sz+(index_sz*2);
	gpu_bytes+=(cat_index_sz*6);
	gpu_bytes+=(comp_sz*2)+src_sz+merge_depth_sz;
	gpu_bytes+=(vert_sz*3);
	std::cout<<"GPU space: "<<gpu_bytes<<" byte(s)\n";

	ggraph->edge_count=edge_count;
	ggraph->vert_count=vert_count;
//	H_ERR(cudaSetDeviceFlags(cudaDeviceMapHost));
	H_ERR(cudaMalloc((void **)&ggraph->adj_card_d, 	index_sz));
	H_ERR(cudaMalloc((void **)&ggraph->beg_pos_d,	index_sz));
	H_ERR(cudaMalloc((void **)&ggraph->adj_list_d, 	edge_sz));
	H_ERR(cudaMalloc((void **)&ggraph->depth_merge, merge_depth_sz));
	H_ERR(cudaMemcpy(ggraph->beg_pos_d,ggraph->beg_pos,index_sz, 
								cudaMemcpyHostToDevice));
	H_ERR(cudaMemcpy(ggraph->adj_card_d,ggraph->adj_card,index_sz, 
								cudaMemcpyHostToDevice));
	H_ERR(cudaMemcpy(ggraph->adj_list_d,ggraph->adj_list,edge_sz,
								cudaMemcpyHostToDevice));
		
	H_ERR(cudaMalloc((void **)&ggraph->ex_sml_off, cat_index_sz));
	H_ERR(cudaMalloc((void **)&ggraph->ex_mid_off, cat_index_sz));
	H_ERR(cudaMalloc((void **)&ggraph->ex_lrg_off, cat_index_sz));
	H_ERR(cudaMalloc((void **)&ggraph->ex_cat_sml_sz,	cat_index_sz));
	H_ERR(cudaMalloc((void **)&ggraph->ex_cat_mid_sz,	cat_index_sz));
	H_ERR(cudaMalloc((void **)&ggraph->ex_cat_lrg_sz,	cat_index_sz));
	
	H_ERR(cudaMalloc((void **)&ggraph->ex_sml_q,	vert_sz));
	H_ERR(cudaMalloc((void **)&ggraph->ex_mid_q,	vert_sz));
	H_ERR(cudaMalloc((void **)&ggraph->ex_lrg_q,	vert_sz));
		
	H_ERR(cudaMalloc((void **)&ggraph->src_list,src_sz));
	H_ERR(cudaMemcpy(ggraph->src_list,src_list,src_sz,cudaMemcpyHostToDevice));
	H_ERR(cudaMalloc((void **)&ggraph->depth_comp_last,comp_sz));
	H_ERR(cudaMalloc((void **)&ggraph->depth_comp_curr,comp_sz));
	cudaStreamCreate(&ggraph->gstream);
	
	for(index_t i=0;i<num_agg_bfs;i++){
		//+----------------------
		//|FOR CLASSIFICATION
		//+----------------------------
	//	H_ERR(cudaMallocHost((void **)&gdata[i]->ex_sml_sz,sizeof(index_t)));
	//	H_ERR(cudaMallocHost((void **)&gdata[i]->ex_mid_sz,sizeof(index_t)));
	//	H_ERR(cudaMallocHost((void **)&gdata[i]->ex_lrg_sz,sizeof(index_t)));
		H_ERR(cudaHostAlloc((void **)&gdata[i]->ex_sml_sz,
					sizeof(index_t),cudaHostAllocMapped));
		H_ERR(cudaHostAlloc((void **)&gdata[i]->ex_mid_sz,
					sizeof(index_t),cudaHostAllocMapped));
		H_ERR(cudaHostAlloc((void **)&gdata[i]->ex_lrg_sz,
					sizeof(index_t),cudaHostAllocMapped));
		H_ERR(cudaHostGetDevicePointer((void **)&gdata[i]->ex_sml_sz_d,
					(void*)gdata[i]->ex_sml_sz,0));
		H_ERR(cudaHostGetDevicePointer((void **)&gdata[i]->ex_mid_sz_d,
					(void*)gdata[i]->ex_mid_sz,0));
		H_ERR(cudaHostGetDevicePointer((void **)&gdata[i]->ex_lrg_sz_d,
					(void*)gdata[i]->ex_lrg_sz,0));
		
		
		gdata[i]->stream=(cudaStream_t *)malloc(sizeof(cudaStream_t)*Q_CARD);
		for(index_t j=0;j<Q_CARD; j++) cudaStreamCreate(&gdata[i]->stream[j]);
	}
	return 0;
}

