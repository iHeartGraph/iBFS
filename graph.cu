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

#include "graph.h"

inline bool deep_than_two(
	vertex_t src,
	depth_t *depth,
	index_t *beg_pos,
	vertex_t *csr
){
	depth[src]=0;
	
	//two level DFS
	index_t beg=beg_pos[src];
	index_t end=beg_pos[src+1];
	for(;beg<end;++beg)
	{
		vertex_t nebr=csr[beg];
		if(depth[nebr] == INFTY)
		{
			depth[nebr]=1;
			index_t nbeg=beg_pos[nebr];
			index_t nend=beg_pos[nebr+1];
			for(;nbeg<nend;++nbeg)
			{
				vertex_t nnebr=csr[nbeg];
				if(depth[nnebr]==INFTY)
					return true;
			}
		}
	}
	return false;
}


graph::graph(
		const char *beg_file,
		const char *csr_file,
		const index_t src_count
):src_count(src_count)
{
	double tm=wtime();
	
	typedef uint64_t index_tt;
	typedef uint64_t vertex_tt;

	//typedef uint32_t index_tt;
	//typedef uint32_t vertex_tt;

	src_list = new vertex_t[src_count];
	vert_count=fsize(beg_file)/sizeof(index_tt) - 1;
	edge_count=fsize(csr_file)/sizeof(vertex_tt);
	FILE *file=fopen(beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<beg_file<<" cannot open\n";
		exit(-1);
	}

	index_tt *tmp_beg_pos = new index_tt[vert_count+1];
	index_tt ret=fread(tmp_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<csr_file<<" cannot open\n";
		exit(-1);
	}

	vertex_tt *tmp_csr = new vertex_tt[edge_count];
	ret=fread(tmp_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);
	
	//converting to uint32_t
	beg_pos = new index_t[vert_count+1];
	csr	= new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		beg_pos[i]=(index_t)tmp_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		csr[i]=(vertex_t)tmp_csr[i];
	
	delete[] tmp_beg_pos;
	delete[] tmp_csr;

	std::cout<<"Graph load (success): "<<vert_count<<" verts, "
		<<edge_count<<" edges "<<wtime()-tm<<" second(s)\n";
}

void
graph::
gen_src()
{
	depth_t *depth=new depth_t[vert_count];
	for(index_t i=0;i<src_count;++i)
	{
		src_list[i] = rand()%vert_count;
		memset(depth, INFTY, sizeof(depth_t)*vert_count);

		//no orphan
		if(!deep_than_two(src_list[i],depth,beg_pos,csr))
		{
			--i;
			continue;
		}

		//no duplication
		for(index_t j=0;j<i;++j)
			if(src_list[i]==src_list[j])
			{
				--i;
				break;
			}
	}

	delete[] depth;
}
