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

//Hang Dec/10/2013
//Hang Mar/05/2015
#include "gpu_ibfs.cuh"
#include "util.h"

__device__ void __sync_warp(int predicate){
	while((!__all(predicate))){
		;
	}
}

//This kernel is executed by one thread
__global__ void init_expand_sort
(	
	vertex_t 	*src_list,
	index_t 	*beg_pos_d,
	index_t		src_grp_off,
	index_t		joint_count,
	index_t		bit_count,
	index_t 	concurr_count,
	comp_t		*depth_comp_last,
	comp_t		*depth_comp_curr,
	index_t		*ex_lrg_sz_d,
	vertex_t	*ex_lrg_q_d
)
{
	index_t tid=threadIdx.x+blockIdx.x*blockDim.x;
	const index_t GRNTY=blockDim.x*gridDim.x;

	while(tid<concurr_count){
		index_t inst_off=tid/bit_count;
		vertex_t vert_id=src_list[tid+src_grp_off];
		index_t	card=beg_pos_d[vert_id+1]-beg_pos_d[vert_id];
		if(!card){return;}
		
		//put all src into big queue 
		//-since we have enough threads 
		ex_lrg_q_d[tid]=vert_id;
	
		comp_t status=(make_uint4(1))<<(tid%bit_count);
		depth_comp_last[vert_id*joint_count+inst_off]=status;
		depth_comp_curr[vert_id*joint_count+inst_off]=status;

		tid+=GRNTY;
	}
	return ;
}

void 
gpu_ibfs::
init_bfs(index_t src_grp_off){
	init_expand_sort
	<<<BLKS_NUM,THDS_NUM>>>
	(
		src_list_d,
		beg_pos_d,
		src_grp_off,
		joint_count,
		bit_count,
		concurr_count,
		depth_comp_last,
		depth_comp_curr,
		ex_lrg_sz_d,
		ex_lrg_q_d
	);
}

//+----------------------
//|for ex_sml_q_d expansion
//+---------------------------
__global__ void td_expand_thd
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	depth_t			curr_level,
	index_t			ex_sml_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz		= ex_sml_sz;
	const index_t GRNLTY	= blockDim.x * gridDim.x;
	index_t tid				= threadIdx.x+blockIdx.x*blockDim.x;

	//used for prefetching
	vertex_t 	ex_ver_curr;
	index_t 	end_curr;
	index_t 	beg_curr;
	vertex_t 	adj_vert;
	comp_t 	vert_depth, adj_depth;

	while(tid<q_sz)
	{
		ex_ver_curr	= ex_q_d[tid];
		beg_curr		= beg_pos_d[ex_ver_curr];
		end_curr		= beg_pos_d[ex_ver_curr+1];
		vert_depth	= depth_comp_last[ex_ver_curr];
		while(beg_curr<end_curr)
		{
			adj_vert=adj_list_d[beg_curr];
			adj_depth=depth_comp_curr[adj_vert];

			if((~(adj_depth.x))&vert_depth.x)
				atomicOr(&(depth_comp_curr[adj_vert].x),vert_depth.x);

			if((~(adj_depth.y))&vert_depth.y)
				atomicOr(&(depth_comp_curr[adj_vert].y),vert_depth.y);

			if((~(adj_depth.z))&vert_depth.z)
				atomicOr(&(depth_comp_curr[adj_vert].z),vert_depth.z);

			if((~(adj_depth.w))&vert_depth.w)
				atomicOr(&(depth_comp_curr[adj_vert].w),vert_depth.w);
			beg_curr++;
		}
		tid += GRNLTY;
	}
}

//+-----------------------
//|ex_mid_q_d expansion
//+------------------------------
__global__ void td_expand_warp
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	depth_t			curr_level,
	index_t			ex_mid_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz		= ex_mid_sz;
	const index_t vec_sz	= 32;
	const index_t TID			= threadIdx.x+blockIdx.x*blockDim.x;
	const index_t lane_s	= TID & (vec_sz-1);
	const index_t GRNLTY	= (blockDim.x*gridDim.x)/vec_sz;
	index_t	vec_id				= TID/vec_sz;

	vertex_t 	ex_ver_curr,adj_vert;
	index_t 	beg_curr, end_curr;
	comp_t 	vert_depth, adj_depth;
	
	while(vec_id<q_sz)
	{
		ex_ver_curr	= ex_q_d[vec_id];
		beg_curr		= beg_pos_d[ex_ver_curr];
		end_curr		= beg_pos_d[ex_ver_curr+1];
		vert_depth	= depth_comp_last[ex_ver_curr];
		index_t lane= beg_curr+lane_s;
		while(lane<end_curr)
		{
			adj_vert=adj_list_d[lane];
			adj_depth=depth_comp_curr[adj_vert];
			if((~(adj_depth.x))&vert_depth.x)
				atomicOr(&(depth_comp_curr[adj_vert].x),vert_depth.x);

			if((~(adj_depth.y))&vert_depth.y)
				atomicOr(&(depth_comp_curr[adj_vert].y),vert_depth.y);

			if((~(adj_depth.z))&vert_depth.z)
				atomicOr(&(depth_comp_curr[adj_vert].z),vert_depth.z);

			if((~(adj_depth.w))&vert_depth.w)
				atomicOr(&(depth_comp_curr[adj_vert].w),vert_depth.w);
			lane+=vec_sz;
		}
		vec_id += GRNLTY;
	}
}

//+-----------------------
//|ex_lrg_q_d expansion
//+------------------------------
__global__ void td_expand_cta
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	depth_t			curr_level,
	index_t			ex_lrg_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz	= ex_lrg_sz;
	index_t	vec_id		= blockIdx.x;

	vertex_t 	ex_ver_curr,adj_vert;
	index_t 	end_curr;
	index_t 	beg_curr;
	comp_t 	vert_depth, adj_depth;
	
	while(vec_id<q_sz)
	{
		ex_ver_curr	= ex_q_d[vec_id];
		beg_curr		= beg_pos_d[ex_ver_curr];
		end_curr		= beg_pos_d[ex_ver_curr+1];
		vert_depth	= depth_comp_last[ex_ver_curr];
		index_t lane= beg_curr+threadIdx.x;
		while(lane<end_curr)
		{
			adj_vert=adj_list_d[lane];
			adj_depth=depth_comp_curr[adj_vert];

			if((~(adj_depth.x))&vert_depth.x)
				atomicOr(&(depth_comp_curr[adj_vert].x),vert_depth.x);

			if((~(adj_depth.y))&vert_depth.y)
				atomicOr(&(depth_comp_curr[adj_vert].y),vert_depth.y);

			if((~(adj_depth.z))&vert_depth.z)
				atomicOr(&(depth_comp_curr[adj_vert].z),vert_depth.z);

			if((~(adj_depth.w))&vert_depth.w)
				atomicOr(&(depth_comp_curr[adj_vert].w),vert_depth.w);
			lane+=blockDim.x;
		}
		vec_id += gridDim.x;
	}
}


//+------------------------
//|ex_mid_q_d expansion
//+------------------------------
__global__ void sw_expand_warp
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	depth_t			curr_level,
	index_t			ex_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz		= ex_sz;
	const index_t vec_sz	= 32;
	const index_t TID			= threadIdx.x+blockIdx.x*blockDim.x;
	const index_t lane_s	= TID & (vec_sz-1);
	const index_t GRNLTY	= (blockDim.x*gridDim.x)/vec_sz;
	index_t	vec_id				= TID/vec_sz;

	vertex_t 	ex_ver_curr;
	index_t 	end_curr, end_curr_revised;
	index_t 	beg_curr;
	unsigned int vert_depth_up, vert_depth_low;
	comp_t 		vert_depth;

	while(vec_id<q_sz)
	{
		ex_ver_curr	= ex_q_d[vec_id];
		beg_curr		= beg_pos_d[ex_ver_curr];
		end_curr		= beg_pos_d[ex_ver_curr+1];
		vert_depth	= depth_comp_curr[ex_ver_curr];
		index_t lane= beg_curr+lane_s;
		
		//make sure all 32 threads can get into the voting loop 
		if((end_curr-beg_curr)&(vec_sz-1))
      end_curr_revised = beg_curr+((((end_curr-beg_curr)>>5)+1)<<5);
    else 
      end_curr_revised = end_curr;

		while(lane<end_curr_revised)
		{
			if(__any((~vert_depth)==0)) break;
			if(lane<end_curr)
				vert_depth|=depth_comp_last[adj_list_d[lane]];
			
			for(int i=16;i>=1;i>>=1){
				vert_depth.x|=__shfl_xor((int)vert_depth.x,i,32);
				vert_depth.y|=__shfl_xor((int)vert_depth.y,i,32);
				vert_depth.z|=__shfl_xor((int)vert_depth.z,i,32);
				vert_depth.w|=__shfl_xor((int)vert_depth.w,i,32);
			}

			lane+=vec_sz;
		}
		
		if(lane_s==0) depth_comp_curr[ex_ver_curr]=vert_depth;
		vec_id += GRNLTY;
	}
}

//+------------------------
//|ex_mid_q_d expansion
//+------------------------------
__global__ void sw_expand_thd
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	depth_t			curr_level,
	index_t			ex_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz		= ex_sz;
	index_t tid			= threadIdx.x+blockIdx.x*blockDim.x;
	const index_t GRNLTY	= blockDim.x*gridDim.x;

	vertex_t 	ex_ver_curr;
	index_t 	beg_curr, end_curr;
	comp_t 	adj_depth_curr, vert_depth;
	
	while(tid<q_sz)
	{
		ex_ver_curr	= ex_q_d[tid];
		end_curr		= beg_pos_d[ex_ver_curr+1];
		beg_curr		= beg_pos_d[ex_ver_curr];
		index_t lane= beg_curr;
		vert_depth=depth_comp_curr[ex_ver_curr];
		while(lane<end_curr)
		{
			if((~vert_depth)==0) break;
			adj_depth_curr=depth_comp_last[adj_list_d[lane]];
			if(((~vert_depth)&adj_depth_curr) != 0) 
				vert_depth|=adj_depth_curr;			
			lane++;
		}

		if((vert_depth^depth_comp_curr[ex_ver_curr]) != 0) 
			depth_comp_curr[ex_ver_curr]=vert_depth;
		tid += GRNLTY;
	}
}

//+-----------------
//|CLFY_EXPAND_SORT
//+----------------------
void 
gpu_ibfs::
td_expand(depth_t	curr_level)
{
	td_expand_thd
	<<<BLKS_NUM, THDS_NUM, 0, stream[0]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_sml_sz[0],
		ex_sml_q_d,
		beg_pos_d,
		adj_list_d
	);

	td_expand_warp
	<<<BLKS_NUM, THDS_NUM, 0, stream[1]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_mid_sz[0],
		ex_mid_q_d,
		beg_pos_d,
		adj_list_d
	);

	
	td_expand_cta
	<<<BLKS_NUM, THDS_NUM, 0, stream[2]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_lrg_sz[0],
		ex_lrg_q_d,
		beg_pos_d,
		adj_list_d
	);

	for(index_t i=0;i<Q_CARD; i++)
		cudaStreamSynchronize(stream[i]);
}

//+----------------------
//|CLFY_EXPAND_SORT
//+----------------------
void 
gpu_ibfs::
bu_expand(depth_t	curr_level)
{
	sw_expand_thd
	<<<BLKS_NUM, THDS_NUM, 0, stream[0]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_sml_sz[0],
		ex_sml_q_d,
		beg_pos_d,
		adj_list_d
	);
		
	sw_expand_thd
	<<<BLKS_NUM, THDS_NUM, 0, stream[1]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_mid_sz[0],
		ex_mid_q_d,
		beg_pos_d,
		adj_list_d
	);
	
	sw_expand_thd
	<<<BLKS_NUM, THDS_NUM, 0, stream[2]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_lrg_sz[0],
		ex_lrg_q_d,
		beg_pos_d,
		adj_list_d
	);

	for(index_t i=0;i<Q_CARD; i++)
		cudaStreamSynchronize(stream[i]);
}

//+----------------------
//|CLFY_EXPAND_SORT
//+----------------------
void 
gpu_ibfs::
sw_expand(depth_t curr_level)
{
	sw_expand_thd
	<<<BLKS_NUM<<1, THDS_NUM, 0, stream[0]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_sml_sz[0],
		ex_sml_q_d,
		beg_pos_d,
		adj_list_d
	);
		
	sw_expand_warp
	<<<BLKS_NUM<<1, THDS_NUM, 0, stream[1]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_mid_sz[0],
		ex_mid_q_d,
		beg_pos_d,
		adj_list_d
	);
	
	sw_expand_warp
	<<<BLKS_NUM<<1, THDS_NUM, 0, stream[2]>>>
	(
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		ex_lrg_sz[0],
		ex_lrg_q_d,
		beg_pos_d,
		adj_list_d
	);

	for(index_t i=0;i<Q_CARD; i++)
		cudaStreamSynchronize(stream[i]);
}

void 
gpu_ibfs::
expander(depth_t level)
{
	if(level<=sw_level)	
		td_expand(level);
	else if((level==sw_level+1))
		sw_expand(level);
	else{//also use merged fq
		//differ sw_expand of no imbalance issue
		bu_expand(level);
	}
}
