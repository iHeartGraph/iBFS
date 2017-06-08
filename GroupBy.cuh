#ifndef __GROUPBY__
#define __GROUPBY__

#include "comm.h"

template<typename vertex_t, typename index_t>
int GroupBy(
	index_t num_bfs,
	vertex_t *source_list,
	vertex_t *new_list,
	csr_graph ggraph
){
	index_t *hubmap=new index_t[ggraph->vert_count];
	index_t *huboff=new index_t[ggraph->vert_count];
	memset(hubmap,0,sizeof(index_t)*ggraph->vert_count);
	memset(huboff,0,sizeof(index_t)*ggraph->vert_count);

	/*counting eligible instances*/
	for(index_t i=0;i<num_bfs;i++)
	{
		index_t beg=ggraph->beg_pos[source_list[i]];
		index_t end=ggraph->beg_pos[source_list[i]+1];
		
		/*out-degree <= sdegree*/
		if(end-beg>128) continue; 
		
		for(index_t j=beg;j<end;j++)
		{
			/*has hub neighbor*/
			if(ggraph->beg_pos[ggraph->adj_list[j]+1]-
				ggraph->beg_pos[ggraph->adj_list[j]] > 128)
			{
				hubmap[ggraph->adj_list[j]]++;
				break;
			}
		}
	}
	
	/*prefix-sum counted data*/
	huboff[0]=0;
	for(index_t i=1;i<ggraph->vert_count;i++)
		huboff[i]=huboff[i-1]+hubmap[i-1];
	index_t remainoff=huboff[ggraph->vert_count-1]+hubmap[ggraph->vert_count-1];
	
	
	index_t test=0;
	/*formulating group*/
	for(index_t i=0;i<num_bfs;i++)
	{
		index_t beg=ggraph->beg_pos[source_list[i]];
		index_t end=ggraph->beg_pos[source_list[i]+1];
		
		/*non-eligible instance*/
		if(end-beg>128)
		{
			new_list[remainoff]=source_list[i];
			remainoff++;
			test++;
			continue;
			
		}
		
		/*eligible instance*/
		bool isEligible=false;
		for(index_t j=beg;j<end;j++)
		{
			/*has hub neighbor*/
			if(ggraph->beg_pos[ggraph->adj_list[j]+1]-
				ggraph->beg_pos[ggraph->adj_list[j]] > 128)
			{
				new_list[huboff[ggraph->adj_list[j]]]=source_list[i];
				huboff[ggraph->adj_list[j]]++;
				isEligible=true;
				break;
			}
		}
		
		if(!isEligible)
		{
			new_list[remainoff]=source_list[i];
			remainoff++;
			test++;
		}
	}
	
	std::ofstream myfile("origin_source.dat");
	std::ofstream newfile("new_source.dat");
	for(index_t i=0;i<num_bfs;i++)
	{
		myfile<<source_list[i]<<"\n";
		newfile<<new_list[i]<<"\n";
	}
	myfile.close();
	newfile.close();


	index_t remain=num_bfs;
	index_t grouped=0;
	std::cout<<"Eligible-BFS(by-source-remain): ";
	for(index_t i=0;i<ggraph->vert_count;i++)
	{
		if(hubmap[i])
		{
			std::cout<<hubmap[i]<<" ";
			grouped+=hubmap[i];
		}
	}
	remain-=grouped;
	std::cout<<remain<<"\n";
	
	std::cout<<"Get test: "<<test<<"\n";
	delete[] hubmap;
	delete[] huboff;
	return 0;
}

#endif
