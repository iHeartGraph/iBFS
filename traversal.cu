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

#include "gpu_ibfs.cuh"

void gpu_ibfs::
debugger(index_t grp_beg)
{
	comp_t *last = new comp_t[g->vert_count*joint_count];
	comp_t *curr = new comp_t[g->vert_count*joint_count];
	index_t *front_count=new index_t[concurr_count];
	
	H_ERR(cudaMemcpy(last, depth_comp_last, 
		sizeof(comp_t)*g->vert_count*joint_count, cudaMemcpyDeviceToHost));
	H_ERR(cudaMemcpy(curr, depth_comp_curr, 
		sizeof(comp_t)*g->vert_count*joint_count, cudaMemcpyDeviceToHost));
	
	//ASSUMING one var can hold all bits
	for(index_t i=0;i<concurr_count;++i)
		front_count[i] = 0;
	
	for(index_t i=0;i<g->vert_count;++i)
	{
		if(last[i]!=curr[i])
		{
			comp_t curr_var=curr[i];
			comp_t last_var=last[i];
			for(index_t j=0;j<concurr_count;++j)
			{
				if((curr_var&(make_uint4(1)<<j))!=(last_var&(make_uint4(1)<<j)))
					++front_count[j];
			}
		}
	}
	
	char filename[256];
	for(index_t i=0;i<concurr_count;++i)
	{
		sprintf(filename, "%d.log",g->src_list[i+grp_beg]);
		FILE *file=fopen(filename, "a");
		fprintf(file, "%d\n",front_count[i]);
		fclose(file);
	}

	delete[] front_count;
	delete[] last;
	delete[] curr;
}


//traversal
void gpu_ibfs::
traversal(index_t	grp_beg){
	cudaSetDevice(0);
	ex_sml_sz[0]	= 0;
	ex_mid_sz[0]	= 0;
	ex_lrg_sz[0]	= concurr_count;
	depth_t level=0;
	
	init_bfs(grp_beg);
	cudaDeviceSynchronize();

	#ifdef ENABLE_MONITORING
	std::cout<<"Vert-ranger: "
		<<g->src_list[grp_beg]<<"~"
		<<g->src_list[grp_beg+concurr_count-1]<<"\n";
	double tm_insp;
	double tm_expd;
	double tm_expand	= 0.0;
	double tm_inspect	= 0.0;
	debugger(grp_beg);
	#endif
	
	is_done[0]=false;	
	for(level=1;;level++){
		if(is_done[0])break;
		
		//----------------------	
		//Expand ex_qs and mark frontier in depth_d	
		//------------------------------------------------
		#ifdef ENABLE_MONITORING
		tm_expd=wtime();	
		#endif
		expander(level);
		cudaDeviceSynchronize();

		#ifdef ENABLE_MONITORING
		tm_expd=wtime()-tm_expd;	
		debugger(grp_beg);
		tm_insp=wtime();	
		#endif
		//-----------------------
		//Generate ex_qs from depth_d
		//---------------------------------
		inspector(level);
		cudaDeviceSynchronize();
		#ifdef ENABLE_MONITORING
		tm_insp=wtime()-tm_insp;
		tm_expand		+= tm_expd;
		tm_inspect	+= tm_insp;
		std::cout<<"@level-"<<(int)level<<"-insp-expd: "
			<<tm_insp<<" "<<tm_expd<<"\n";
		#endif
	}
	
	#ifdef ENABLE_MONITORING
	debugger(grp_beg);
	std::cout<<"insp-expd: "<<tm_inspect<<" "<<tm_expand<<"\n";
	#endif
}
