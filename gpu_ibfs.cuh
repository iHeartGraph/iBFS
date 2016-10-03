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

//Graph format: Json based format: [src_id, src_weigh,[[connected_ver_0, edge_weight],[connected_ver_1, edge_weight],[connected_ver_2, edge_weight]]]
//Storage format: 
//struct{
//		int: src_ver
//		Arr: [ver_0|ver_1|ver_2|...]
//		Int: num_conn_ver
//	}
#ifndef	GPU_IBFS_H
#define	GPU_IBFS_H

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include "util.h"
#include "wtime.h"
#include "graph.h"

class gpu_ibfs{
	
	//variable
public:
	const graph *g;
	
	vertex_t* src_list_d;
	vertex_t*	adj_list_d;
	index_t*	beg_pos_d;
	vertex_t*	ex_sml_q_d;
	vertex_t*	ex_mid_q_d;
	vertex_t*	ex_lrg_q_d;
	index_t*	cat_sml_sz_d;
	index_t*	cat_mid_sz_d;
	index_t*	cat_lrg_sz_d;
	index_t*	cat_sml_off_d;
	index_t*	cat_mid_off_d;
	index_t*	cat_lrg_off_d;
	
	index_t	*ex_sml_sz;
	index_t	*ex_mid_sz;
	index_t	*ex_lrg_sz;
	index_t	*ex_sml_sz_d;
	index_t	*ex_mid_sz_d;
	index_t	*ex_lrg_sz_d;
	
	comp_t	*depth_comp_last;
	comp_t	*depth_comp_curr;
	cudaStream_t *stream, gstream;
	
	index_t	sml_shed;
	index_t	lrg_shed;
	index_t	gpu_count;
	index_t	*gpu_ranger;
	index_t joint_count;
	index_t bit_count;
	index_t concurr_count;
	index_t	sw_level;
	bool *is_done;
	bool *is_done_d;

	//constructor
public:
	gpu_ibfs(){};
	gpu_ibfs(	
		const graph *g,
		index_t concurr_count,
		depth_t			sw_level,
		index_t			gpu_count,
		index_t			sml_shed,
		index_t			lrg_shed);
	~gpu_ibfs(){};

	//functions
public:
	int write_result();
	
	void traversal(index_t src_grp_off);
	void init_bfs(index_t src_grp_off);
	void debugger(index_t);

	//expander series
	void expander(depth_t level);
	void bu_expand(depth_t level);
	void td_expand(depth_t level);
	void sw_expand(depth_t level);

	//inspector series
	void inspector(depth_t level);
	void td_inspect(depth_t level);
	void sw_inspect(depth_t level);
};
#endif
