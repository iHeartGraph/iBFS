//Dec/15/2013
//Feb/16/2015
//Mar/4/2015
#include "comm.h"


template<typename depth_t, typename comp_t, typename index_t>
__global__ void extract_depth(
	comp_t *depth_comp_last,
	comp_t *depth_comp_curr,
	depth_t *depth_merge,
	index_t num_agg_bfs,
	depth_t	curr_level,
	index_t	vert_count
){
	
	//+-
	//|Diff of comp_last and comp_curr = just visited verts
	//+--------------------------------------
	index_t tid=threadIdx.x+blockIdx.x*blockDim.x;
	index_t wid=tid>>5;
	const index_t vec_sz=32;
	const index_t LANE=tid%32;
	const index_t GRNTY=blockDim.x*gridDim.x;
	const index_t WGRNTY=GRNTY>>5;
	const depth_t LEVEL=curr_level;
	const index_t COMP_COUNT=vert_count;
	comp_t comp_last, comp_curr;
//	unsigned int recv_up, recv_low;
	unsigned int change[4], recv;
		
	while(tid<COMP_COUNT){
		comp_last=depth_comp_last[tid];
		comp_curr=depth_comp_curr[tid];
		change[0]=comp_curr.x-comp_last.x;
		change[1]=comp_curr.y-comp_last.y;
		change[2]=comp_curr.z-comp_last.z;
		change[3]=comp_curr.w-comp_last.w;
		
		if(change[0] || change[1] || change[2] || change[3]){
			depth_comp_last[tid]=comp_curr;

			for(int i=0;i<vec_sz;i++){
				for(int j=0;j<4;j++){
					recv=change[j];
					recv=__shfl((int)recv,i);
					if(recv&(1<<LANE)) 
						depth_merge[wid*vec_sz*num_agg_bfs+i*num_agg_bfs+j*vec_sz+LANE]=LEVEL;
				}
			}
		}
		__syncthreads();
		tid += GRNTY;
		wid += WGRNTY;
	}
}


template
<	typename vertex_t,
	typename index_t,
	typename depth_t>
__global__ void td_frontier_count 
(
	index_t	 	*ex_cat_sml_sz,
	index_t	 	*ex_cat_mid_sz,
	index_t	 	*ex_cat_lrg_sz,
	comp_t 		*depth_comp_last,
	comp_t 		*depth_comp_curr,
	depth_t 	*depth_merge,
	index_t		*adj_card_d,
	index_t		num_agg_bfs,
	depth_t	 	curr_level,
	index_t	 	vert_count,
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	//+-
	//|Diff of comp_last and comp_curr = just visited verts
	//+--------------------------------------
	index_t TID_ST=threadIdx.x+blockIdx.x*blockDim.x;
	index_t tid=TID_ST;
	index_t wid=TID_ST>>5;
	const index_t vec_sz=32;
	const index_t LANE=tid%32;
	const index_t GRNTY=blockDim.x*gridDim.x;
	const index_t WGRNTY=GRNTY>>5;
	const depth_t LEVEL=curr_level;
	const index_t COMP_COUNT=vert_count;
	comp_t comp_last, comp_curr;
	index_t card_curr;
	unsigned int recv;
	
	index_t ex_sml_ct= 0;
	index_t ex_mid_ct= 0;
	index_t ex_lrg_ct= 0;
	unsigned int change[4];

	while(tid<COMP_COUNT){
		comp_last=depth_comp_last[tid];
		comp_curr=depth_comp_curr[tid];
		change[0]=comp_curr.x-comp_last.x;
		change[1]=comp_curr.y-comp_last.y;
		change[2]=comp_curr.z-comp_last.z;
		change[3]=comp_curr.w-comp_last.w;
		
		if(change[0] || change[1] || change[2] || change[3]){
			
			//depth_comp_last[tid]=comp_curr;
			//CANNOT CHANGE NOW!!!
			//NEED diff FOR GATHERING
			for(int i=0;i<vec_sz;i++){
				for(int j=0;j<4;j++){
					recv=change[j];
					recv=__shfl((int)recv,i);
					if(recv&(1<<LANE)) 
						depth_merge[wid*vec_sz*num_agg_bfs+i*num_agg_bfs+j*vec_sz+LANE]=LEVEL;
				}
			}

			card_curr=adj_card_d[tid];
			if(card_curr>0)
				if(card_curr<sml_shed) ex_sml_ct++;
				else if(card_curr>lrg_shed) ex_lrg_ct++;
				else ex_mid_ct++;
		}

		__syncthreads();
		tid += GRNTY;
		wid += WGRNTY;
	}
	
	ex_cat_sml_sz[TID_ST]= ex_sml_ct;
 	ex_cat_mid_sz[TID_ST]= ex_mid_ct;
 	ex_cat_lrg_sz[TID_ST]= ex_lrg_ct;
}

template <typename vertex_t, typename index_t,
			typename depth_t>
__global__ void sw_frontier_count 
(
	index_t	 *ex_cat_sml_sz,
	index_t	 *ex_cat_mid_sz,
	index_t	 *ex_cat_lrg_sz,
	comp_t 	*depth_comp_curr,
	index_t		*adj_card_d,
	index_t		num_agg_bfs,
	depth_t	 curr_level,
	index_t	 vert_count,
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	const index_t TID_ST	= threadIdx.x +blockIdx.x*blockDim.x;
	const index_t NUM_VER	= vert_count;	
	const index_t GRNTY		= blockDim.x*gridDim.x;

	index_t card_curr,card_next;
	comp_t depth_curr,depth_next;
	index_t ex_sml_ct	= 0;
	index_t ex_mid_ct	= 0;
	index_t ex_lrg_ct	= 0;
	index_t tid				= TID_ST;
	
	//Figure out its own start and end pos
	//We want each thread to inspect a continuous block
	index_t step_sz = NUM_VER/GRNTY;
	if(step_sz<16){
		//small problem size
		const index_t REMAINDER = NUM_VER-step_sz*GRNTY;
		if(TID_ST<REMAINDER) step_sz ++;
		const index_t beg_pos=step_sz*TID_ST+ 
					(TID_ST >= REMAINDER ? REMAINDER:0);
		const index_t end_pos=beg_pos+step_sz;
		tid=beg_pos;
		if(step_sz){
			card_curr		= adj_card_d[beg_pos];
			depth_curr	= depth_comp_curr[beg_pos];
//			depth_comp_last[beg_pos]=depth_curr;
		}

		while(tid<end_pos){
			tid++;
			if(tid<end_pos){
				card_next		= adj_card_d[tid];
				depth_next	= depth_comp_curr[tid];
//				depth_comp_last[tid]=depth_next;
			}
				
			if((card_curr>0)&&(
				depth_curr.x!=0xffffffff ||
				depth_curr.y!=0xffffffff ||
				depth_curr.z!=0xffffffff ||
				depth_curr.w!=0xffffffff 
				)){
				if(card_curr<sml_shed) ex_sml_ct++;
				else if(card_curr>lrg_shed) ex_lrg_ct++;
				else ex_mid_ct++;
			}
			
			depth_curr=depth_next;
			card_curr	= card_next;
		}
	}else{
		//big problem size
		//we want each thread to get 16x indices to check.
		if(step_sz&0xf) step_sz=((step_sz>>4)+1)<<4;
		if(NUM_VER-step_sz*GRNTY>0) step_sz<<=1;	
		__shared__ index_t	beg_pos_s[THDS_NUM];
		__shared__ depth_t	mark_s[THDS_NUM<<4];
		index_t	beg_pos_own			= step_sz*TID_ST;
		beg_pos_s[threadIdx.x] 	= beg_pos_own;
		const index_t TRIES 		= step_sz>>4;
		const index_t lane_id		= threadIdx.x%32;
		const index_t	warp_id		= threadIdx.x>>5;
		const index_t	thread_off= threadIdx.x<<4;
		index_t tries = 0;
		index_t	load_ptr;
		index_t	proc_vert;

		//+--------------------
		//in shared memory data representation
		//---------------------------------------
		//			0000|1111
		//					 ||||
		//					 |||+-------- INFTY ? 
		//					 ||+--------- SML   ?
		//					 |+---------- MID   ?
		//					 +----------- LRG   ?
		//FURTHER OPTIMIZATION CAN BE CARRIED OUT
		//BY PACKING MORE HERE.
		//-----------------------------------------
		while(tries<TRIES){
			//warp stride loading
			for(int i=0; i<32; i+=2){
				proc_vert=(warp_id<<5)+i+(lane_id>>4);
				depth_curr.x=0xffffffff;
				depth_curr.y=0xffffffff;
				depth_curr.z=0xffffffff;
				depth_curr.w=0xffffffff;
				if(proc_vert<NUM_VER)
				{
					load_ptr=(tries<<4)+beg_pos_s[proc_vert]+(lane_id%16);
					if(load_ptr<NUM_VER)
					{
						depth_curr=depth_comp_curr[load_ptr];
						card_curr=adj_card_d[load_ptr];
					//	depth_comp_last[load_ptr]=depth_curr;
					}else card_curr=0;//out-of-boundary, set it orphan
				}else card_curr=0;//out-of-boundary, set it orphan
				mark_s[(proc_vert<<4)+(lane_id%16)]=0x00;

				//construct hub-cache
//				if(depth_curr == LEVEL){
//					hub_ptr	= load_ptr & (HUB_BU_SZ-1);
//					if(card_curr > hub_card[hub_ptr]){
//						hub_vert[hub_ptr] = load_ptr;
//						hub_card[hub_ptr] = card_curr;
//					}
//				}else
				if((card_curr>0) &&(
				depth_curr.x!=0xffffffff ||
				depth_curr.y!=0xffffffff ||
				depth_curr.z!=0xffffffff ||
				depth_curr.w!=0xffffffff))
					if(card_curr<sml_shed)
						mark_s[(proc_vert<<4)+(lane_id%16)]=0x03;//0011
					else if(card_curr>lrg_shed)
						mark_s[(proc_vert<<4)+(lane_id%16)]=0x09;//1001
					else mark_s[(proc_vert<<4)+(lane_id%16)]=0x05;//0101
			}
			__syncthreads();

			//thread stride checking
			for(int i=0; i<16; i++)
				if(mark_s[thread_off+i]&0x02){
					ex_sml_ct++;
				}else if(mark_s[thread_off+i]&0x08){
					ex_lrg_ct++;
				}else if(mark_s[thread_off+i]&0x04){
					ex_mid_ct++;
				}
			++tries;
			__syncthreads();
		}
	}
	
	ex_cat_sml_sz[TID_ST]= ex_sml_ct;
	ex_cat_mid_sz[TID_ST]= ex_mid_ct;
	ex_cat_lrg_sz[TID_ST]= ex_lrg_ct;
}

template <typename vertex_t, typename index_t,
			typename depth_t>
__global__ void td_frontier_gather 
(
	index_t	 	*ex_sml_q,
	index_t	 	*ex_mid_q,
	index_t	 	*ex_lrg_q,
	index_t	 	*ex_sml_off,
	index_t	 	*ex_mid_off,
	index_t	 	*ex_lrg_off,
	comp_t 		*depth_comp_last,
	comp_t 		*depth_comp_curr,
	depth_t 	*depth_merge,
	index_t		*adj_card_d,
	index_t		num_agg_bfs,
	depth_t	 	curr_level,
	index_t	 	vert_count,
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	index_t TID_ST=threadIdx.x+blockIdx.x*blockDim.x;
	index_t tid=TID_ST;
	const index_t GRNLTY=blockDim.x*gridDim.x;
	const index_t COMP_COUNT=vert_count;
	comp_t comp_last, comp_curr;
	index_t card_curr;
	
	index_t ex_sml_ptr	= ex_sml_off[TID_ST];
	index_t ex_mid_ptr	= ex_mid_off[TID_ST];
	index_t ex_lrg_ptr	= ex_lrg_off[TID_ST];
	unsigned int change[4];

	
	while(tid<COMP_COUNT){
		comp_last=depth_comp_last[tid];
		comp_curr=depth_comp_curr[tid];
		change[0]=comp_curr.x-comp_last.x;
		change[1]=comp_curr.y-comp_last.y;
		change[2]=comp_curr.z-comp_last.z;
		change[3]=comp_curr.w-comp_last.w;
		
		if(change[0] || change[1] || change[2] || change[3]){		
//			depth_comp_last[tid].x=comp_curr.x;
//			depth_comp_last[tid].y=comp_curr.y;
//			depth_comp_last[tid].z=comp_curr.z;
//			depth_comp_last[tid].w=comp_curr.w;
			depth_comp_last[tid]=comp_curr;
			card_curr=adj_card_d[tid];
			if(card_curr>0)
				if(card_curr<sml_shed)
				{
					ex_sml_q[ex_sml_ptr]=tid;
					ex_sml_ptr++;
				}
				else if(card_curr>lrg_shed)
				{
					ex_lrg_q[ex_lrg_ptr]=tid;
					ex_lrg_ptr++;
				}
				else
				{
					ex_mid_q[ex_mid_ptr]=tid;
					ex_mid_ptr++;
				}
		}
		tid+= GRNLTY;
	}
}


template <typename vertex_t, typename index_t,
			typename depth_t>
__global__ void sw_frontier_gather 
(
	index_t	 	*ex_sml_q,
	index_t	 	*ex_mid_q,
	index_t	 	*ex_lrg_q,
	index_t	 	*ex_sml_off,
	index_t	 	*ex_mid_off,
	index_t	 	*ex_lrg_off,
	comp_t 	*depth_comp_curr,
	index_t		*adj_card_d,
	index_t		num_agg_bfs,
	depth_t	 curr_level,
	index_t	 vert_count,
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	const index_t TID_ST	= threadIdx.x +blockIdx.x*blockDim.x;
	const index_t NUM_VER	= vert_count;	
	const index_t GRNTY		= blockDim.x*gridDim.x;

	index_t ex_sml_ptr	= ex_sml_off[TID_ST];
	index_t ex_mid_ptr	= ex_mid_off[TID_ST];
	index_t ex_lrg_ptr	= ex_lrg_off[TID_ST];

	index_t card_curr, card_next;
	comp_t depth_curr,depth_next;
	index_t	tid,tid_next;
	
	//Figure out its own start and end pos
	//We want each thread to inspect a continuous block
	index_t step_sz = NUM_VER/GRNTY;
	if(step_sz<16){
		//small problem size
		const index_t REMAINDER = NUM_VER-step_sz*GRNTY;
		if(TID_ST<REMAINDER) step_sz ++;
		const index_t beg_pos=step_sz*TID_ST+ 
					(TID_ST >= REMAINDER ? REMAINDER:0);
		const index_t end_pos=beg_pos+step_sz;
		
		tid=beg_pos;tid_next=beg_pos;
		if(step_sz){
			card_curr		= adj_card_d[beg_pos];
			depth_curr	= depth_comp_curr[beg_pos];
			tid_next++;
		}
			
		while(tid<end_pos){
			if(tid_next < end_pos){
				card_next		= adj_card_d[tid_next];
				depth_next	= depth_comp_curr[tid_next];
			}
			
			if((card_curr>0)&&(
				depth_curr.x!=0xffffffff ||
				depth_curr.y!=0xffffffff ||
				depth_curr.z!=0xffffffff ||
				depth_curr.w!=0xffffffff 
				)){
				if(card_curr<sml_shed)
				{
					ex_sml_q[ex_sml_ptr]=tid;
					ex_sml_ptr++;
				}
				else if(card_curr>lrg_shed)
				{
					ex_lrg_q[ex_lrg_ptr]=tid;
					ex_lrg_ptr++;
				}
				else
				{
					ex_mid_q[ex_mid_ptr]=tid;
					ex_mid_ptr++;
				}
			}
			tid=tid_next;
			tid_next++;
			card_curr	= card_next;
			depth_curr=depth_next;
		}
	}else{
		//big problem size
		//We want each thread to get 16x indices to check.
		if(step_sz&0xf) step_sz=((step_sz>>4)+1)<<4;
		if(NUM_VER-step_sz*GRNTY>0) step_sz<<=1;	
		__shared__ index_t	beg_pos_s[THDS_NUM];
		__shared__ depth_t	mark_s[THDS_NUM<<4];
		index_t	beg_pos_own			= step_sz*TID_ST;
		beg_pos_s[threadIdx.x] 	= beg_pos_own;
		const index_t TRIES 		= step_sz>>4;
		const index_t lane_id		= threadIdx.x%32;
		const index_t	warp_id		= threadIdx.x>>5;
		const index_t	thread_off= threadIdx.x<<4;
		index_t tries = 0;
		index_t	load_ptr;
		index_t	proc_vert;

		//+--------------------
		//in shared memory data representation
		//---------------------------------------
		//			0000|1111
		//					 ||||
		//					 |||+-------- INFTY ? 
		//					 ||+--------- SML   ?
		//					 |+---------- MID   ?
		//					 +----------- LRG   ?
		//FURTHER OPTIMIZATION CAN BE CARRIED OUT
		//BY PACKING MORE HERE.
		//-----------------------------------------
		while(tries<TRIES){
			//warp stride loading
			for(int i=0; i<32; i+=2){
				proc_vert=(warp_id<<5)+i+(lane_id>>4);
				depth_curr.x=0xffffffff;
				depth_curr.y=0xffffffff;
				depth_curr.z=0xffffffff;
				depth_curr.w=0xffffffff;
				if(proc_vert<NUM_VER)
				{
					load_ptr=(tries<<4)+beg_pos_s[proc_vert]+(lane_id%16);
					if(load_ptr<NUM_VER)
					{
						depth_curr=depth_comp_curr[load_ptr];
						card_curr=adj_card_d[load_ptr];
					}else card_curr=0;//out-of-boundary, set it orphan
				}else card_curr=0;//out-of-boundary, set it orphan
				mark_s[(proc_vert<<4)+(lane_id%16)]=0x00;

				if((card_curr>0) &&(
				depth_curr.x!=0xffffffff ||
				depth_curr.y!=0xffffffff ||
				depth_curr.z!=0xffffffff ||
				depth_curr.w!=0xffffffff))
					if(card_curr<sml_shed)
						mark_s[(proc_vert<<4)+(lane_id%16)]=0x03;//0011
					else if(card_curr>lrg_shed)
						mark_s[(proc_vert<<4)+(lane_id%16)]=0x09;//1001
					else mark_s[(proc_vert<<4)+(lane_id%16)]=0x05;//0101
			}
			__syncthreads();

			//thread stride checking
			for(int i=0; i<16; i++)
				if(mark_s[thread_off+i]&0x02){
					ex_sml_q[ex_sml_ptr]=beg_pos_own+(tries<<4)+i;
					ex_sml_ptr++;
				}else if(mark_s[thread_off+i]&0x08){
					ex_lrg_q[ex_lrg_ptr]=beg_pos_own+(tries<<4)+i;
					ex_lrg_ptr++;
				}else if(mark_s[thread_off+i]&0x04){
					ex_mid_q[ex_mid_ptr]=beg_pos_own+(tries<<4)+i;
					ex_mid_ptr++;
				}

			++tries;
			__syncthreads();
		}
	}
}


//+--------------------------------------------
//|CLASSIFIED BASED STORAGE FOR EX_QUEUE
//+---------------------------------------------------
template< typename vertex_t, 
			typename index_t,
			typename depth_t>
void td_inspect
(
	tdata *gdata,
	csr_graph ggraph,
	index_t		num_agg_bfs,
	index_t	 sw_level,
	depth_t	 curr_level,
	index_t  vert_count,//num ver in the graph
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	td_frontier_count	
	<vertex_t, index_t, depth_t>
	<<<BLKS_NUM,THDS_NUM>>>
	(
		ggraph->ex_cat_sml_sz,
		ggraph->ex_cat_mid_sz,
		ggraph->ex_cat_lrg_sz,
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		ggraph->depth_merge,
		ggraph->adj_card_d,
		num_agg_bfs,
		curr_level,
		vert_count,
		sml_shed,
		lrg_shed
	);
	cudaDeviceSynchronize();
	
	insp_scan
 	<vertex_t, index_t>
 	(
 		ggraph->ex_cat_sml_sz,
 		ggraph->ex_sml_off,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		gdata[0]->ex_sml_sz_d,
 		gdata[0]->stream[0]
 	);
 	insp_scan
 	<vertex_t, index_t>
 	(
 		ggraph->ex_cat_mid_sz,
 		ggraph->ex_mid_off,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		gdata[0]->ex_mid_sz_d,
 		gdata[0]->stream[1]
 	);
 	insp_scan
 	<vertex_t, index_t>
 	(
 		ggraph->ex_cat_lrg_sz,
 		ggraph->ex_lrg_off,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		gdata[0]->ex_lrg_sz_d,
 		gdata[0]->stream[2]
 	);
 	

	cudaStreamSynchronize(gdata[0]->stream[0]);
	cudaStreamSynchronize(gdata[0]->stream[1]);
	cudaStreamSynchronize(gdata[0]->stream[2]);
	cudaDeviceSynchronize();
	
	td_frontier_gather	
	<vertex_t, index_t, depth_t>
	<<<BLKS_NUM,THDS_NUM>>>
	(
		ggraph->ex_sml_q,
		ggraph->ex_mid_q,
		ggraph->ex_lrg_q,
		ggraph->ex_sml_off,
		ggraph->ex_mid_off,
		ggraph->ex_lrg_off,
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		ggraph->depth_merge,
		ggraph->adj_card_d,
		num_agg_bfs,
		curr_level,
		vert_count,
		sml_shed,
		lrg_shed
	);

}

template< typename vertex_t, 
			typename index_t,
			typename depth_t>
void sw_inspect
(
	tdata *gdata,
	csr_graph ggraph,
	index_t		num_agg_bfs,
	depth_t	 curr_level,
	index_t  vert_count,//num ver in the graph
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	ENABLE_BTUP=true;
	extract_depth
	<depth_t, comp_t, index_t>
	<<<BLKS_NUM<<2,THDS_NUM<<2,0,ggraph->gstream>>>
	(
		ggraph->depth_comp_last,
		ggraph->depth_comp_curr,
		ggraph->depth_merge,
		num_agg_bfs,
		curr_level,
		vert_count
	);
	
	//ONLY COMPRESS ONE TIME
//	if(ENABLE_BTUP==false){
//		ENABLE_BTUP=true;
//		comp_depth
//		<depth_t, comp_t, index_t>
//		<<<BLKS_NUM<<2,THDS_NUM<<2>>>
//		(
//			ggraph->depth_merge,
//			ggraph->depth_comp_last,
//			ggraph->depth_comp_curr,
//			curr_level,
//			0,
//			vert_count
//		);
//		cudaDeviceSynchronize();
//	}
	
	sw_frontier_count
 	<vertex_t,index_t,depth_t>
 	<<<BLKS_NUM,THDS_NUM>>>
 	(
 		ggraph->ex_cat_sml_sz,
 		ggraph->ex_cat_mid_sz,
 		ggraph->ex_cat_lrg_sz,
 		ggraph->depth_comp_curr,
 		ggraph->adj_card_d,
 		num_agg_bfs,
 		curr_level,
 		vert_count,
 		sml_shed,
 		lrg_shed
 	);
	
	cudaDeviceSynchronize();
	
	insp_scan
 	<vertex_t, index_t>
 	(
 		ggraph->ex_cat_sml_sz,
 		ggraph->ex_sml_off,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		gdata[0]->ex_sml_sz_d,
 		gdata[0]->stream[0]
 	);
 	insp_scan
 	<vertex_t, index_t>
 	(
 		ggraph->ex_cat_mid_sz,
 		ggraph->ex_mid_off,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		gdata[0]->ex_mid_sz_d,
 		gdata[0]->stream[1]
 	);
 	insp_scan
 	<vertex_t, index_t>
 	(
 		ggraph->ex_cat_lrg_sz,
 		ggraph->ex_lrg_off,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		gdata[0]->ex_lrg_sz_d,
 		gdata[0]->stream[2]
 	);
 	

	cudaStreamSynchronize(gdata[0]->stream[0]);
	cudaStreamSynchronize(gdata[0]->stream[1]);
	cudaStreamSynchronize(gdata[0]->stream[2]);
	cudaDeviceSynchronize();
	
	sw_frontier_gather
	<vertex_t, index_t, depth_t>
	<<<BLKS_NUM, THDS_NUM>>>
	(
		ggraph->ex_sml_q,
		ggraph->ex_mid_q,
		ggraph->ex_lrg_q,
		ggraph->ex_sml_off,
		ggraph->ex_mid_off,
		ggraph->ex_lrg_off,
 		ggraph->depth_comp_curr,
		ggraph->adj_card_d,
		num_agg_bfs,
		curr_level,
		vert_count,
		sml_shed,
		lrg_shed
	);
	cudaStreamSynchronize(ggraph->gstream);
}


template<typename vertex_t, 
			typename index_t,
			typename depth_t,
			typename comp_t>
void inspector
(
	tdata *gdata,
	csr_graph ggraph,
	index_t		num_agg_bfs,
	index_t		sw_level,
	depth_t	 	level,
	index_t  	vert_count,
	const index_t sml_shed,
	const index_t lrg_shed
)
{
//		comp_t *d_status_curr,*d_status_last;
//		cudaMallocHost((void **)&d_status_curr, sizeof(comp_t)*vert_count);
//		cudaMallocHost((void **)&d_status_last, sizeof(comp_t)*vert_count);
//		cudaMemcpy(d_status_curr, ggraph->depth_comp_curr,sizeof(comp_t)*vert_count,cudaMemcpyDeviceToHost);
//		cudaMemcpy(d_status_last, ggraph->depth_comp_last,sizeof(comp_t)*vert_count,cudaMemcpyDeviceToHost);
//		
//		int ft_count=0; int diff=0;
//		for(int i=0;i<vert_count;i++){
//			if(d_status_curr[i]!=d_status_last[i]) diff++;
//			if(d_status_curr[i]!=0xffffffffffffffff && ggraph->adj_card[i]>0) ft_count++;
//		}
//		std::cout<<"ft_count-diff "<<ft_count<<" "<<diff<<"\n";
//		cudaFreeHost(d_status_last);
//		cudaFreeHost(d_status_curr);
//	
	
	if(level<sw_level)
		td_inspect
		<vertex_t, index_t, depth_t>
		(
			gdata,
			ggraph,
			num_agg_bfs,
			sw_level,
			level,
			vert_count,
			sml_shed,
			lrg_shed
		);
	else{
		
		sw_inspect
		<vertex_t, index_t, depth_t>
		(
			gdata,
			ggraph,
			num_agg_bfs,
			level,
			vert_count,
			sml_shed,
			lrg_shed
		);
	}
}
