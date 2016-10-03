/*
 * Copyright 2016 The George Washington University
 * Written by Hang Liu 
 * Directed by Prof. Howie Huang
 *
 * https://www.seas.gwu.edu/~howie/
 * Contact: iheartgraph@gmail.com
 *
 * 
 * Please cite the following paper:
 * 
 * Hang Liu, H. Howie Huang and Yang Hu. 2016. iBFS: Concurrent Breadth-First Search on GPUs. Proceedings of the 2016 International Conference on Management of Data. ACM. 
 *
 * This file is part of iBFS.
 *
 * iBFS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * iBFS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with iBFS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "scan.cuh"
#include "wtime.h"
#include "validate.h"
#include "gpu_ibfs.cuh"

//constructor
gpu_ibfs::
gpu_ibfs(
	const graph* g,
	index_t concurr_count,
	depth_t sw_level,
	index_t gpu_count,
	index_t sml_shed,
	index_t lrg_shed)
	:g(g),sw_level(sw_level),gpu_count(gpu_count),
	sml_shed(sml_shed),lrg_shed(lrg_shed),concurr_count(concurr_count)
{
	bit_count = sizeof(comp_t)<<3;
	joint_count = (index_t)ceil((concurr_count*1.0)/bit_count);
	
	std::cout<<"joint_count vs bit_count: "<<joint_count<<" "<<bit_count<<"\n";
	cudaSetDevice(0);	
	long gpu_bytes= 0;
	
	const size_t cat_index_sz 	= sizeof(index_t)*BLKS_NUM*THDS_NUM;
	const size_t edge_sz 				= sizeof(vertex_t)*g->edge_count;
	const size_t vert_sz 				= sizeof(vertex_t)*g->vert_count;
	const	size_t index_sz 			= sizeof(index_t)*g->vert_count;
	const	size_t comp_sz 				= sizeof(comp_t)*g->vert_count*joint_count;
	const size_t src_sz					= sizeof(vertex_t)*g->src_count;

	H_ERR(cudaMalloc((void **)&beg_pos_d,		index_sz));
	H_ERR(cudaMalloc((void **)&adj_list_d, 	edge_sz));
	gpu_bytes	+= edge_sz+(index_sz);
		
	H_ERR(cudaMalloc((void **)&cat_sml_off_d, cat_index_sz));
	H_ERR(cudaMalloc((void **)&cat_mid_off_d, cat_index_sz));
	H_ERR(cudaMalloc((void **)&cat_lrg_off_d, cat_index_sz));
	H_ERR(cudaMalloc((void **)&cat_sml_sz_d,	cat_index_sz));
	H_ERR(cudaMalloc((void **)&cat_mid_sz_d,	cat_index_sz));
	H_ERR(cudaMalloc((void **)&cat_lrg_sz_d,	cat_index_sz));
	gpu_bytes+=(cat_index_sz*6);
	
	H_ERR(cudaMalloc((void **)&ex_sml_q_d,	vert_sz));
	H_ERR(cudaMalloc((void **)&ex_mid_q_d,	vert_sz));
	H_ERR(cudaMalloc((void **)&ex_lrg_q_d,	vert_sz));
	gpu_bytes+=(vert_sz*3);
		
	H_ERR(cudaMalloc((void **)&src_list_d,src_sz));
	H_ERR(cudaMalloc((void **)&depth_comp_last,comp_sz));
	H_ERR(cudaMalloc((void **)&depth_comp_curr,comp_sz));
	cudaStreamCreate(&gstream);
	
	gpu_bytes+=(comp_sz*2)+src_sz;
	
	H_ERR(cudaHostAlloc((void **)&is_done,
			sizeof(bool),cudaHostAllocMapped));
	H_ERR(cudaHostGetDevicePointer((void **)&is_done_d,
			(void*)is_done,0));
	
	//+----------------------
	//|FOR CLASSIFICATION
	//+----------------------------
	H_ERR(cudaHostAlloc((void **)&ex_sml_sz,
				sizeof(index_t),cudaHostAllocMapped));
	H_ERR(cudaHostAlloc((void **)&ex_mid_sz,
				sizeof(index_t),cudaHostAllocMapped));
	H_ERR(cudaHostAlloc((void **)&ex_lrg_sz,
				sizeof(index_t),cudaHostAllocMapped));
	H_ERR(cudaHostGetDevicePointer((void **)&ex_sml_sz_d,
				(void*)ex_sml_sz,0));
	H_ERR(cudaHostGetDevicePointer((void **)&ex_mid_sz_d,
				(void*)ex_mid_sz,0));
	H_ERR(cudaHostGetDevicePointer((void **)&ex_lrg_sz_d,
				(void*)ex_lrg_sz,0));
	
	stream=(cudaStream_t *)malloc(sizeof(cudaStream_t)*Q_CARD);
	for(index_t j=0;j<Q_CARD; j++) cudaStreamCreate(&stream[j]);
	std::cout<<"GPU space: "<<gpu_bytes<<" byte(s)\n";

	//copy graph data + source list to GPU
	H_ERR(cudaMemcpy(beg_pos_d,g->beg_pos,index_sz, 
			cudaMemcpyHostToDevice));
	H_ERR(cudaMemcpy(adj_list_d,g->csr,edge_sz,
			cudaMemcpyHostToDevice));
	H_ERR(cudaMemcpy(src_list_d,g->src_list,src_sz,
			cudaMemcpyHostToDevice));
}

