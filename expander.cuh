//Hang Dec/10/2013
//Hang Mar/05/2015

#include "comm.h"

__device__ void __sync_warp(int predicate){
	while((!__all(predicate))){
		;
	}
}

//This kernel is executed by one thread
template<typename vertex_t, 
		typename index_t,
		typename depth_t>
__global__ void init_expand_sort
(	
	vertex_t 	*src_list,
	index_t		src_grp_off,
	depth_t		*depth_merge,
	comp_t		*depth_comp_last,
	comp_t		*depth_comp_curr,
	index_t		num_agg_bfs,
	vertex_t	*parent_d,
	index_t		*adj_card_d,
	index_t		*ex_sml_sz_d,
	index_t		*ex_mid_sz_d,
	index_t		*ex_lrg_sz_d,
	vertex_t	*ex_sml_q,
	vertex_t	*ex_mid_q,
	vertex_t	*ex_lrg_q,
	index_t		vert_count,
	index_t		sml_shed,
	index_t		lrg_shed
)
{
//	parent_d[src_v]	=-1;
	if(threadIdx.x==num_agg_bfs){
		ex_sml_sz_d[0]	= 0;
		ex_mid_sz_d[0]	= 0;
		ex_lrg_sz_d[0]	= num_agg_bfs;
	}

	if(threadIdx.x<num_agg_bfs){
		vertex_t src_v=src_list[num_agg_bfs*src_grp_off+threadIdx.x];
	//	depth_merge[src_v*num_agg_bfs+threadIdx.x]= 0;
	//	
	//	index_t	card	= adj_card_d[src_v];
	//	if(!card) return;
		
		//put all src into big queue since we have enough threads 
		//initally
		ex_lrg_q[threadIdx.x]=src_v;
		
		unsigned int status=1;
		if(threadIdx.x<32){
			status<<=threadIdx.x;
			depth_comp_last[src_v].x=status;
			depth_comp_curr[src_v].x=status;
		}else if(threadIdx.x<64){
			status<<=(threadIdx.x-32);
			depth_comp_last[src_v].y=status;
			depth_comp_curr[src_v].y=status;
		}else if(threadIdx.x<96){
			status<<=(threadIdx.x-64);
			depth_comp_last[src_v].z=status;
			depth_comp_curr[src_v].z=status;
		}else if(threadIdx.x<128){
			status<<=(threadIdx.x-96);
			depth_comp_last[src_v].w=status;
			depth_comp_curr[src_v].w=status;
		}
	}
	
	return ;
}

template<typename vertex_t,
				typename index_t,
				typename depth_t>
void init_bfs(
	vertex_t	*src_v, 
	index_t		src_grp_off,
	tdata 		*gdata,
	csr_graph ggraph,
	index_t 	num_agg_bfs,
	index_t		sml_shed,
	index_t		lrg_shed
){

	init_expand_sort
	<vertex_t, index_t, depth_t>
	<<<1,THDS_NUM>>>
	(
		ggraph->src_list,
		src_grp_off,
		ggraph->depth_merge,
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		ggraph->adj_card_d,
		gdata[0]->ex_sml_sz_d,
		gdata[0]->ex_mid_sz_d,
		gdata[0]->ex_lrg_sz_d,
		ggraph->ex_sml_q,
		ggraph->ex_mid_q,
		ggraph->ex_lrg_q,
		ggraph->vert_count,
		sml_shed,
		lrg_shed
	);

//cudaDeviceSynchronize();
//	index_t vert_count=ggraph->vert_count;
//		comp_t *d_status_curr,*d_status_last;
//		cudaMallocHost((void **)&d_status_curr, sizeof(comp_t)*vert_count);
//		cudaMallocHost((void **)&d_status_last, sizeof(comp_t)*vert_count);
//		cudaMemcpy(d_status_curr, ggraph->depth_comp_curr,sizeof(comp_t)*vert_count,cudaMemcpyDeviceToHost);
//		cudaMemcpy(d_status_last, ggraph->depth_comp_last,sizeof(comp_t)*vert_count,cudaMemcpyDeviceToHost);
//		
//		int ft_count=0; int diff=0;
//		for(int i=0;i<vert_count;i++){
//			if(d_status_curr[i]!=0) diff++;
//			if(d_status_curr[i]!=0xffffffffffffffff && ggraph->adj_card[i]>0) ft_count++;
//		}
//		std::cout<<"ft_count-diff "<<ft_count<<" "<<diff<<"\n";
//		cudaFreeHost(d_status_last);
//		cudaFreeHost(d_status_curr);
//	
}

//+----------------------
//|for ex_sml_q expansion
//+---------------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t,
		typename comp_t>
__global__ void td_expand_thd
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	index_t			num_agg_bfs,
	vertex_t		*parent_d,
	depth_t			curr_level,
	index_t			ex_sml_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*	__restrict__ adj_card_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz		= ex_sml_sz;
	const index_t GRNLTY	= blockDim.x * gridDim.x;
	index_t tid				= threadIdx.x+blockIdx.x*blockDim.x;

	//used for prefetching
	vertex_t 	ex_ver_curr;
	index_t 	card_curr;
	index_t 	beg_pos_curr;
	vertex_t 	adj_vert;
	comp_t vert_depth,adj_depth;
	unsigned int vert_depth_vec[4];

	while(tid<q_sz)
	{
		ex_ver_curr	= ex_q_d[tid];
		card_curr		= adj_card_d[ex_ver_curr];
		beg_pos_curr= beg_pos_d[ex_ver_curr];
		index_t lane= beg_pos_curr;
		card_curr+= beg_pos_curr;
		vert_depth=depth_comp_last[ex_ver_curr];
		vert_depth_vec[0]=vert_depth.x;
		vert_depth_vec[1]=vert_depth.y;
		vert_depth_vec[2]=vert_depth.z;
		vert_depth_vec[3]=vert_depth.w;
		
		while(lane<card_curr)
		{
			adj_vert=adj_list_d[lane];
			adj_depth=depth_comp_curr[adj_vert];

			if((~(adj_depth.x))&vert_depth_vec[0])
				atomicOr(&(depth_comp_curr[adj_vert].x),vert_depth_vec[0]);

			if((~(adj_depth.y))&vert_depth_vec[1])
				atomicOr(&(depth_comp_curr[adj_vert].y),vert_depth_vec[1]);

			if((~(adj_depth.z))&vert_depth_vec[2])
				atomicOr(&(depth_comp_curr[adj_vert].z),vert_depth_vec[2]);

			if((~(adj_depth.w))&vert_depth_vec[3])
				atomicOr(&(depth_comp_curr[adj_vert].w),vert_depth_vec[3]);
			lane++;
		}
		tid += GRNLTY;
	}
}

//+-----------------------
//|ex_mid_q expansion
//+------------------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t,
		typename comp_t>
__global__ void td_expand_warp
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	index_t			num_agg_bfs,
	vertex_t		*parent_d,
	depth_t			curr_level,
	index_t			ex_mid_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*	__restrict__ adj_card_d,
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
	index_t 	card_curr;
	index_t 	beg_pos_curr;
	comp_t 	vert_depth,adj_depth;
	unsigned int vert_depth_vec[4];
	
	while(vec_id<q_sz)
	{
		ex_ver_curr	= ex_q_d[vec_id];
		beg_pos_curr= beg_pos_d[ex_ver_curr];
		card_curr		= adj_card_d[ex_ver_curr];
		vert_depth=depth_comp_last[ex_ver_curr];
		vert_depth_vec[0]=vert_depth.x;
		vert_depth_vec[1]=vert_depth.y;
		vert_depth_vec[2]=vert_depth.z;
		vert_depth_vec[3]=vert_depth.w;
		index_t lane= beg_pos_curr+lane_s;
		card_curr+= beg_pos_curr;
		while(lane<card_curr)
		{
			adj_vert=adj_list_d[lane];
			adj_depth=depth_comp_curr[adj_vert];

			if((~(adj_depth.x))&vert_depth_vec[0])
				atomicOr(&(depth_comp_curr[adj_vert].x),vert_depth_vec[0]);

			if((~(adj_depth.y))&vert_depth_vec[1])
				atomicOr(&(depth_comp_curr[adj_vert].y),vert_depth_vec[1]);

			if((~(adj_depth.z))&vert_depth_vec[2])
				atomicOr(&(depth_comp_curr[adj_vert].z),vert_depth_vec[2]);

			if((~(adj_depth.w))&vert_depth_vec[3])
				atomicOr(&(depth_comp_curr[adj_vert].w),vert_depth_vec[3]);
			lane+=vec_sz;
		}
		vec_id += GRNLTY;
	}
}

//+-----------------------
//|ex_lrg_q expansion
//+------------------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t,
		typename comp_t>
__global__ void td_expand_cta
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	index_t			num_agg_bfs,
	vertex_t		*parent_d,
	depth_t			curr_level,
	index_t			ex_lrg_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*	__restrict__ adj_card_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz	= ex_lrg_sz;
	index_t	vec_id		= blockIdx.x;

	vertex_t 	ex_ver_curr,adj_vert;
	index_t 	card_curr;
	index_t 	beg_pos_curr;
	comp_t 	vert_depth,adj_depth;
	unsigned int vert_depth_vec[4];
	
	while(vec_id<q_sz)
	{
		ex_ver_curr	= ex_q_d[vec_id];
		beg_pos_curr= beg_pos_d[ex_ver_curr];
		card_curr		= adj_card_d[ex_ver_curr];
		vert_depth=depth_comp_last[ex_ver_curr];
		vert_depth_vec[0]=vert_depth.x;
		vert_depth_vec[1]=vert_depth.y;
		vert_depth_vec[2]=vert_depth.z;
		vert_depth_vec[3]=vert_depth.w;
		index_t lane= beg_pos_curr+threadIdx.x;
		card_curr+= beg_pos_curr;
		while(lane<card_curr)
		{
			adj_vert=adj_list_d[lane];
			adj_depth=depth_comp_curr[adj_vert];

			if((~(adj_depth.x))&vert_depth_vec[0])
				atomicOr(&(depth_comp_curr[adj_vert].x),vert_depth_vec[0]);

			if((~(adj_depth.y))&vert_depth_vec[1])
				atomicOr(&(depth_comp_curr[adj_vert].y),vert_depth_vec[1]);

			if((~(adj_depth.z))&vert_depth_vec[2])
				atomicOr(&(depth_comp_curr[adj_vert].z),vert_depth_vec[2]);

			if((~(adj_depth.w))&vert_depth_vec[3])
				atomicOr(&(depth_comp_curr[adj_vert].w),vert_depth_vec[3]);
			lane+=blockDim.x;
		}
		vec_id += gridDim.x;
	}
}


//+------------------------
//|ex_mid_q expansion
//+------------------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t,
		typename comp_t>
__global__ void sw_expand_warp
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	index_t			num_agg_bfs,
	vertex_t		*parent_d,
	depth_t			curr_level,
	index_t			ex_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*	__restrict__ adj_card_d,
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
	index_t 	card_curr, card_curr_revised;
	index_t 	beg_pos_curr;
	unsigned int vert_depth_vec[4];
	comp_t 		vert_depth;

	while(vec_id<q_sz)
	{
		ex_ver_curr	= ex_q_d[vec_id];
		beg_pos_curr= beg_pos_d[ex_ver_curr];
		card_curr		= adj_card_d[ex_ver_curr];
		vert_depth=depth_comp_curr[ex_ver_curr];
		vert_depth_vec[0]=vert_depth.x;
		vert_depth_vec[1]=vert_depth.y;
		vert_depth_vec[2]=vert_depth.z;
		vert_depth_vec[3]=vert_depth.w;
		index_t lane= beg_pos_curr+lane_s;
		
		//make sure all 32 threads can get into the voting loop 
		if(card_curr%vec_sz)
      card_curr_revised = (((card_curr>>5)+1)<<5);
    else 
      card_curr_revised = card_curr;

		card_curr+= beg_pos_curr;
		card_curr_revised+=beg_pos_curr;
		while(lane<card_curr_revised)
		{
			/*if(	vert_depth_vec[0]==0xffffffff &&
					vert_depth_vec[1]==0xffffffff &&
					vert_depth_vec[2]==0xffffffff &&
					vert_depth_vec[3]==0xffffffff
					) break;*/
			if(lane<card_curr){
				vert_depth=depth_comp_last[adj_list_d[lane]];
				
				vert_depth_vec[0]|=vert_depth.x;
				vert_depth_vec[1]|=vert_depth.y;
				vert_depth_vec[2]|=vert_depth.z;
				vert_depth_vec[3]|=vert_depth.w;
			}

			for(int i=16;i>=1;i>>=1){
				vert_depth_vec[0]|=__shfl_xor((int)vert_depth_vec[0],i,32);
				vert_depth_vec[1]|=__shfl_xor((int)vert_depth_vec[1],i,32);
				vert_depth_vec[2]|=__shfl_xor((int)vert_depth_vec[2],i,32);
				vert_depth_vec[3]|=__shfl_xor((int)vert_depth_vec[3],i,32);
			}
			lane+=vec_sz;
		}
		
		if(lane_s==0){
			depth_comp_curr[ex_ver_curr].x=vert_depth_vec[0];
			depth_comp_curr[ex_ver_curr].y=vert_depth_vec[1];
			depth_comp_curr[ex_ver_curr].z=vert_depth_vec[2];
			depth_comp_curr[ex_ver_curr].w=vert_depth_vec[3];
		}
		
		vec_id += GRNLTY;
	}
}

//+------------------------
//|ex_mid_q expansion
//+------------------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t,
		typename comp_t>
__global__ void sw_expand_thd
(	
	comp_t			*depth_comp_last,
	comp_t			*depth_comp_curr,
	index_t			num_agg_bfs,
	vertex_t		*parent_d,
	depth_t			curr_level,
	index_t			ex_sz,
	const vertex_t* __restrict__ ex_q_d,
	const index_t*	__restrict__ adj_card_d,
	const index_t*  __restrict__ beg_pos_d,
	const vertex_t* __restrict__ adj_list_d
)
{
	const index_t q_sz		= ex_sz;
	index_t tid			= threadIdx.x+blockIdx.x*blockDim.x;
	const index_t GRNLTY	= blockDim.x*gridDim.x;

	vertex_t 	ex_ver_curr;
	index_t 	card_curr;
	index_t 	beg_pos_curr;
	comp_t 		 vert_depth;
	unsigned int vert_depth_vec[4];
	
	while(tid<q_sz)
	{
		ex_ver_curr	= ex_q_d[tid];
		card_curr		= adj_card_d[ex_ver_curr];
		beg_pos_curr= beg_pos_d[ex_ver_curr];
		index_t lane= beg_pos_curr;
		card_curr+= beg_pos_curr;
		vert_depth=depth_comp_curr[ex_ver_curr];
		vert_depth_vec[0]=vert_depth.x;
		vert_depth_vec[1]=vert_depth.y;
		vert_depth_vec[2]=vert_depth.z;
		vert_depth_vec[3]=vert_depth.w;
		while(lane<card_curr)
		{
			/*if(	vert_depth_vec[0]==0xffffffff &&
					vert_depth_vec[1]==0xffffffff &&
					vert_depth_vec[2]==0xffffffff &&
					vert_depth_vec[3]==0xffffffff
					) break;*/
			vert_depth=depth_comp_last[adj_list_d[lane]];
			
			vert_depth_vec[0]|=vert_depth.x;
			vert_depth_vec[1]|=vert_depth.y;
			vert_depth_vec[2]|=vert_depth.z;
			vert_depth_vec[3]|=vert_depth.w;
			lane++;
		}
		
		depth_comp_curr[ex_ver_curr].x=vert_depth_vec[0];
		depth_comp_curr[ex_ver_curr].y=vert_depth_vec[1];
		depth_comp_curr[ex_ver_curr].z=vert_depth_vec[2];
		depth_comp_curr[ex_ver_curr].w=vert_depth_vec[3];
		tid += GRNLTY;
	}
}

//+-----------------
//|CLFY_EXPAND_SORT
//+----------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t>
void td_expand
(
	tdata	*gdata,
	csr_graph ggraph,
	index_t	num_agg_bfs,
	depth_t	curr_level
)
{
	td_expand_thd<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[0]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_sml_sz[0],
		ggraph->ex_sml_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);

	td_expand_warp<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[1]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_mid_sz[0],
		ggraph->ex_mid_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);

	
	td_expand_cta<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[2]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_lrg_sz[0],
		ggraph->ex_lrg_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);

	for(index_t i=0;i<Q_CARD; i++)
		cudaStreamSynchronize(gdata[0]->stream[i]);
}

//+----------------------
//|CLFY_EXPAND_SORT
//+----------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t>
void bu_expand
(	
	tdata	*gdata,
	csr_graph ggraph,
	index_t	num_agg_bfs,
	depth_t	curr_level
)
{
	sw_expand_thd<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[0]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_sml_sz[0],
		ggraph->ex_sml_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);
		
	sw_expand_thd<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[1]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_mid_sz[0],
		ggraph->ex_mid_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);
	
	sw_expand_thd<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[2]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_lrg_sz[0],
		ggraph->ex_lrg_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);

	for(index_t i=0;i<Q_CARD; i++)
		cudaStreamSynchronize(gdata[0]->stream[i]);
}

//+----------------------
//|CLFY_EXPAND_SORT
//+----------------------
template<typename vertex_t, 
		typename index_t,
		typename depth_t,
		typename comp_t>
void sw_expand
(	
	tdata	*gdata,
	csr_graph ggraph,
	index_t	num_agg_bfs,
	depth_t	curr_level
)
{
	sw_expand_thd<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[0]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_sml_sz[0],
		ggraph->ex_sml_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);
		
	sw_expand_warp<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[1]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_mid_sz[0],
		ggraph->ex_mid_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);
	
	sw_expand_warp<vertex_t,index_t,depth_t,comp_t>
	<<<BLKS_NUM, THDS_NUM, 0, gdata[0]->stream[2]>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		num_agg_bfs,
		gdata[0]->parent_d,
		curr_level,
		gdata[0]->ex_lrg_sz[0],
		ggraph->ex_lrg_q,
		ggraph->adj_card_d,
		ggraph->beg_pos_d,
		ggraph->adj_list_d
	);

	for(index_t i=0;i<Q_CARD; i++)
		cudaStreamSynchronize(gdata[0]->stream[i]);
}


template<typename vertex_t, 
		typename index_t,
		typename depth_t,
		typename comp_t>
void expander(	
	tdata	*gdata,
	csr_graph ggraph,
	index_t	num_agg_bfs,
	index_t	sw_level,
	depth_t	level
)
{
	if(level<=sw_level)	
		td_expand
		<vertex_t, index_t, depth_t>
		(
			gdata,
			ggraph,
			num_agg_bfs,
			level
		);
	else if((level==sw_level+1||level==sw_level+2||level==sw_level+3))
		sw_expand
		<vertex_t, index_t, depth_t,comp_t>
		(
			gdata,
			ggraph,
			num_agg_bfs,
			level
		);
	else{//also use merged fq
		bu_expand//but do not consider imbalance issue 
		<vertex_t, index_t, depth_t>
		(
			gdata,
			ggraph,
			num_agg_bfs,
			level
		);
	}
}
