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

//Dec/15/2013
//Feb/16/2015
//Mar/4/2015
#include "util.h"
#include "gpu_ibfs.cuh"
#include "scan.cuh"

__global__ void extract_depth(
	comp_t *depth_comp_last,
	comp_t *depth_comp_curr,
	bool 			*is_done_d,
	depth_t	curr_level,
	index_t	vert_count
){
	
	//+-
	//|Diff of comp_last and comp_curr = just visited verts
	//+--------------------------------------

	index_t tid=threadIdx.x+blockIdx.x*blockDim.x;
	index_t wid=tid>>5;
//	const index_t vec_sz=32;
//	const index_t LANE=tid%32;
	const index_t GRNTY=blockDim.x*gridDim.x;
	const index_t WGRNTY=GRNTY>>5;
//	const depth_t LEVEL=curr_level;
	const index_t COMP_COUNT=vert_count;
	comp_t comp_last, comp_curr;
//	unsigned int recv_up, recv_low;
	bool is_done=true;

	while(tid<COMP_COUNT){
		comp_last=depth_comp_last[tid];
		comp_curr=depth_comp_curr[tid];
		if(comp_last != comp_curr){
			depth_comp_last[tid]=comp_curr;
			if(is_done) is_done=false;
		}

	//	for(int i=0;i<vec_sz;i++){
	//		recv_up=change>>32;
	//		recv_low=change&0xffffffff;
	//		recv_up=__shfl((int)recv_up,i);
	//		recv_low=__shfl((int)recv_low,i);

	//		if(recv_up&(1<<LANE)) 
	//			depth_merge[wid*vec_sz*num_agg_bfs+i*num_agg_bfs+vec_sz+LANE]=LEVEL;
	//		if(recv_low&(1<<LANE)) 
	//			depth_merge[wid*vec_sz*num_agg_bfs+i*num_agg_bfs+LANE]=LEVEL;
	//	}
	//	__syncthreads();
		tid += GRNTY;
		wid += WGRNTY;
	}
//is_done_d[0]=false;
	if(is_done==false) is_done_d[0]=false;
}

__global__ void td_frontier_count(
	index_t		*beg_pos_d,
	index_t	 	*cat_sml_sz_d,
	index_t	 	*cat_mid_sz_d,
	index_t	 	*cat_lrg_sz_d,
	comp_t 		*depth_comp_last,
	comp_t 		*depth_comp_curr,
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
//	const index_t vec_sz=32;
//	const index_t LANE=tid%32;
	const index_t GRNTY=blockDim.x*gridDim.x;
	const index_t WGRNTY=GRNTY>>5;
//	const depth_t LEVEL=curr_level;
	const index_t COMP_COUNT=vert_count;
	comp_t comp_last, comp_curr;
	index_t card_curr;
//	unsigned int recv_up, recv_low;
	
	index_t ex_sml_ct= 0;
	index_t ex_mid_ct= 0;
	index_t ex_lrg_ct= 0;

	while(tid<COMP_COUNT){
		comp_last=depth_comp_last[tid];
		comp_curr=depth_comp_curr[tid];
		if((comp_last^comp_curr) != 0){
			
			//depth_comp_last[tid]=comp_curr;
			//CANNOT CHANGE NOW!!!
			//NEED diff FOR GATHERING

			card_curr=beg_pos_d[tid+1]-beg_pos_d[tid];
			if(card_curr>0)
				if(card_curr<sml_shed) ex_sml_ct++;
				else if(card_curr>lrg_shed) ex_lrg_ct++;
				else ex_mid_ct++;
		}

	//	for(int i=0;i<vec_sz;i++){
	//		recv_up=change>>32;
	//		recv_low=change&0xffffffff;
	//		recv_up=__shfl((int)recv_up,i);
	//		recv_low=__shfl((int)recv_low,i);

	//		if(recv_up&(1<<LANE)) 
	//			depth_merge[wid*vec_sz*num_agg_bfs+i*num_agg_bfs+vec_sz+LANE]=LEVEL;
	//		if(recv_low&(1<<LANE)) 
	//			depth_merge[wid*vec_sz*num_agg_bfs+i*num_agg_bfs+LANE]=LEVEL;
	//	}
		__syncthreads();
		tid += GRNTY;
		wid += WGRNTY;
	}
	
	cat_sml_sz_d[TID_ST]= ex_sml_ct;
 	cat_mid_sz_d[TID_ST]= ex_mid_ct;
 	cat_lrg_sz_d[TID_ST]= ex_lrg_ct;
}

__global__ void sw_frontier_count(
	index_t		*beg_pos_d,
	index_t	 *cat_sml_sz_d,
	index_t	 *cat_mid_sz_d,
	index_t	 *cat_lrg_sz_d,
	comp_t 	*depth_comp_curr,
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
			card_curr		= beg_pos_d[beg_pos+1]-beg_pos_d[beg_pos];
			depth_curr	= depth_comp_curr[beg_pos];
//			depth_comp_last[beg_pos]=depth_curr;
		}

		while(tid<end_pos){
			tid++;
			if(tid<end_pos){
				card_next		= beg_pos_d[tid+1]-beg_pos_d[tid];
				depth_next	= depth_comp_curr[tid];
//				depth_comp_last[tid]=depth_next;
			}
				
			if(((~depth_curr)!=0) && (card_curr>0)){
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
				depth_curr.x=0;depth_curr.y=0;depth_curr.y=0;depth_curr.w=0;
				depth_curr=~depth_curr;
				if(proc_vert<NUM_VER)
				{
					load_ptr=(tries<<4)+beg_pos_s[proc_vert]+(lane_id%16);
					if(load_ptr<NUM_VER)
					{
						depth_curr=depth_comp_curr[load_ptr];
						card_curr=beg_pos_d[load_ptr+1]-beg_pos_d[load_ptr];
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
				if(((~depth_curr)!=0) && (card_curr>0))
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
	
	cat_sml_sz_d[TID_ST]= ex_sml_ct;
	cat_mid_sz_d[TID_ST]= ex_mid_ct;
	cat_lrg_sz_d[TID_ST]= ex_lrg_ct;
}

__global__ void td_frontier_gather(
	index_t		*beg_pos_d,
	index_t	 	*ex_sml_q_d,
	index_t	 	*ex_mid_q_d,
	index_t	 	*ex_lrg_q_d,
	index_t	 	*cat_sml_off_d,
	index_t	 	*cat_mid_off_d,
	index_t	 	*cat_lrg_off_d,
	comp_t 		*depth_comp_last,
	comp_t 		*depth_comp_curr,
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
	
	index_t ex_sml_ptr	= cat_sml_off_d[TID_ST];
	index_t ex_mid_ptr	= cat_mid_off_d[TID_ST];
	index_t ex_lrg_ptr	= cat_lrg_off_d[TID_ST];
	
	while(tid<COMP_COUNT){
		comp_last=depth_comp_last[tid];
		comp_curr=depth_comp_curr[tid];
		comp_t change=comp_curr-comp_last;
		if(change!=0){
			depth_comp_last[tid]=comp_curr;
			card_curr=beg_pos_d[tid+1]-beg_pos_d[tid];
			if(card_curr>0)
				if(card_curr<sml_shed)
				{
					ex_sml_q_d[ex_sml_ptr]=tid;
					ex_sml_ptr++;
				}
				else if(card_curr>lrg_shed)
				{
					ex_lrg_q_d[ex_lrg_ptr]=tid;
					ex_lrg_ptr++;
				}
				else
				{
					ex_mid_q_d[ex_mid_ptr]=tid;
					ex_mid_ptr++;
				}
		}
		tid+= GRNLTY;
	}
}


__global__ void sw_frontier_gather(
	index_t		*beg_pos_d,
	index_t	 	*ex_sml_q_d,
	index_t	 	*ex_mid_q_d,
	index_t	 	*ex_lrg_q_d,
	index_t	 	*cat_sml_off_d,
	index_t	 	*cat_mid_off_d,
	index_t	 	*cat_lrg_off_d,
	comp_t 	*depth_comp_curr,
	depth_t	 curr_level,
	index_t	 vert_count,
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	const index_t TID_ST	= threadIdx.x +blockIdx.x*blockDim.x;
	const index_t NUM_VER	= vert_count;	
	const index_t GRNTY		= blockDim.x*gridDim.x;

	index_t ex_sml_ptr	= cat_sml_off_d[TID_ST];
	index_t ex_mid_ptr	= cat_mid_off_d[TID_ST];
	index_t ex_lrg_ptr	= cat_lrg_off_d[TID_ST];

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
			card_curr		= beg_pos_d[beg_pos+1]-beg_pos_d[beg_pos];
			depth_curr	= depth_comp_curr[beg_pos];
			tid_next++;
		}
			
		while(tid<end_pos){
			if(tid_next < end_pos){
				card_next		= beg_pos_d[tid_next+1]-beg_pos_d[tid_next];
				depth_next	= depth_comp_curr[tid_next];
			}
			
			if((card_curr>0)&&((~depth_curr)!=0)){
				if(card_curr<sml_shed)
				{
					ex_sml_q_d[ex_sml_ptr]=tid;
					ex_sml_ptr++;
				}
				else if(card_curr>lrg_shed)
				{
					ex_lrg_q_d[ex_lrg_ptr]=tid;
					ex_lrg_ptr++;
				}
				else
				{
					ex_mid_q_d[ex_mid_ptr]=tid;
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
				depth_curr.x=0;depth_curr.y=0;depth_curr.y=0;depth_curr.w=0;
				depth_curr=~depth_curr;
				if(proc_vert<NUM_VER)
				{
					load_ptr=(tries<<4)+beg_pos_s[proc_vert]+(lane_id%16);
					if(load_ptr<NUM_VER)
					{
						depth_curr=depth_comp_curr[load_ptr];
						card_curr=beg_pos_d[load_ptr+1]-beg_pos_d[load_ptr];
					}else card_curr=0;//out-of-boundary, set it orphan
				}else card_curr=0;//out-of-boundary, set it orphan
				mark_s[(proc_vert<<4)+(lane_id%16)]=0x00;

				if(((~depth_curr)!=0)&&(card_curr>0))
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
					ex_sml_q_d[ex_sml_ptr]=beg_pos_own+(tries<<4)+i;
					ex_sml_ptr++;
				}else if(mark_s[thread_off+i]&0x08){
					ex_lrg_q_d[ex_lrg_ptr]=beg_pos_own+(tries<<4)+i;
					ex_lrg_ptr++;
				}else if(mark_s[thread_off+i]&0x04){
					ex_mid_q_d[ex_mid_ptr]=beg_pos_own+(tries<<4)+i;
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
void 
gpu_ibfs::
td_inspect(depth_t curr_level)
{
	td_frontier_count	
	<<<BLKS_NUM,THDS_NUM>>>(
		beg_pos_d,
		cat_sml_sz_d,
		cat_mid_sz_d,
		cat_lrg_sz_d,
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		g->vert_count,
		sml_shed,
		lrg_shed
	);
	cudaDeviceSynchronize();
	
	insp_scan
 	(
 		cat_sml_sz_d,
 		cat_sml_off_d,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		ex_sml_sz_d,
 		stream[0]
 	);
 	insp_scan
 	(
 		cat_mid_sz_d,
 		cat_mid_off_d,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		ex_mid_sz_d,
 		stream[1]
 	);
 	insp_scan
 	(
 		cat_lrg_sz_d,
 		cat_lrg_off_d,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		ex_lrg_sz_d,
 		stream[2]
 	);
 	

	cudaStreamSynchronize(stream[0]);
	cudaStreamSynchronize(stream[1]);
	cudaStreamSynchronize(stream[2]);
	cudaDeviceSynchronize();
	
	td_frontier_gather	
	<<<BLKS_NUM,THDS_NUM>>>(
		beg_pos_d,
		ex_sml_q_d,
		ex_mid_q_d,
		ex_lrg_q_d,
		cat_sml_off_d,
		cat_mid_off_d,
		cat_lrg_off_d,
		depth_comp_last,
		depth_comp_curr,
		curr_level,
		g->vert_count,
		sml_shed,
		lrg_shed
	);
}

void 
gpu_ibfs::
sw_inspect(depth_t curr_level)
{
	is_done[0]=true;
	extract_depth
	<<<BLKS_NUM<<2,THDS_NUM<<2,0,gstream>>>
	(
		depth_comp_last,
		depth_comp_curr,
		is_done_d,
		curr_level,
		g->vert_count
	);
	
	sw_frontier_count
 	<<<BLKS_NUM,THDS_NUM>>>(
		beg_pos_d,
 		cat_sml_sz_d,
 		cat_mid_sz_d,
 		cat_lrg_sz_d,
 		depth_comp_curr,
 		curr_level,
 		g->vert_count,
 		sml_shed,
 		lrg_shed
	);
	
	cudaDeviceSynchronize();
	
	insp_scan
 	(
 		cat_sml_sz_d,
 		cat_sml_off_d,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		ex_sml_sz_d,
 		stream[0]
 	);
 	insp_scan
 	(
 		cat_mid_sz_d,
 		cat_mid_off_d,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		ex_mid_sz_d,
 		stream[1]
 	);
 	insp_scan
 	(
 		cat_lrg_sz_d,
 		cat_lrg_off_d,
 		THDS_NUM*BLKS_NUM,
 		BLKS_NUM>>1,
 		THDS_NUM>>1,
 		ex_lrg_sz_d,
 		stream[2]
 	);

	cudaStreamSynchronize(stream[0]);
	cudaStreamSynchronize(stream[1]);
	cudaStreamSynchronize(stream[2]);
	cudaDeviceSynchronize();
	
	sw_frontier_gather
	<<<BLKS_NUM, THDS_NUM>>>(
		beg_pos_d,
		ex_sml_q_d,
		ex_mid_q_d,
		ex_lrg_q_d,
		cat_sml_off_d,
		cat_mid_off_d,
		cat_lrg_off_d,
 		depth_comp_curr,
		curr_level,
		g->vert_count,
		sml_shed,
		lrg_shed
	);
	cudaStreamSynchronize(gstream);
}

void gpu_ibfs::
inspector(depth_t level)
{
	if(level<sw_level)
		td_inspect(level);
	else 
		sw_inspect(level);
}
