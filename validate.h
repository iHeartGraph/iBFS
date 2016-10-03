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

#include "util.h"
#include <iostream>

template < 	typename index_t,
			typename vertex_t,
			typename depth_t> 
index_t validate(	
		depth_t 	*depth_h, 
		vertex_t	*adj_list,
		index_t		*adj_card,
		index_t		*beg_pos,
		index_t		vert_count)
{
	
	//divide the adjacency vertices of this vertex into two parts
	//-parent part: all the vertices connected to the parent part 
	//				should have smaller or equal to base_depth
	//-children part: all vertices connected to the children part
	//				should have larger or equal to base_depth
	//depth difference should be no more than 1
	depth_t base_depth;
	char err_type;
	for(int i = 0; i< vert_count; i++){
		base_depth = depth_h[i];
		if(base_depth ==INFTY) continue;

		if(!base_depth){//src vertex
			for(int j = beg_pos[i]; j< beg_pos[i]+adj_card[i]; j++)
				if((depth_h[adj_list[j]]!= 1) && (adj_list[j]!=i)){
					std::cout<<"Vert: "<<adj_list[j]
								<<" depth: "<<(int)depth_h[adj_list[j]]<<"\n";
					return -1;		
				}
		}else{
			
			//neighbor should consist both parents and children
			err_type = 0x00;//0x00000011
			//												||
			//												|+----child bit
			//												+-----parent bit
			//
			for(int j = beg_pos[i]; j<beg_pos[i]+adj_card[i]; j++){
				if(depth_h[adj_list[j]] == base_depth - 1)
					if(!(err_type & 0x02)) err_type |=0x02;
				else if(depth_h[adj_list[j]] == base_depth + 1)
					if(!(err_type & 0x01)) err_type |=0x01;
				else if(depth_h[adj_list[j]] != base_depth)
					return -2;
			}
			
			if(err_type == 0x01){
				std::cout<<"Vert: "<<i<<" err: "<<(int)err_type<<"\n";
				return -3;
			}	
		}
	}
	return 0;
}


template< 	typename vertex_t,
			typename index_t,
			typename depth_t>
void report(	
				index_t &agg_tr_edges,
				index_t &agg_tr_v,
				depth_t	*depth_h,
				index_t	*adj_card,
				index_t	vert_count)
{
	agg_tr_edges	= 0;
	agg_tr_v		= 0;
	for(index_t i = 0; i< vert_count; i++)
		if(depth_h[i] != INFTY){
			agg_tr_v ++;
			agg_tr_edges += adj_card[i];
		}
}
