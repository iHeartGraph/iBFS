//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...


//#include "graph.h"
#include <fstream>
#include <omp.h>
#include "GroupBy.cuh"

#define FILE_NOT_EXIST		1
#define FILE_EXIST				0

template<typename index_t>
void progress(
		index_t &report,
		index_t beg_time,
		int			tid,
		index_t prc_line
){
	if(prc_line>report && tid==0){
			std::cout<<prc_line<<" lines\t"
					<<"time: "<<wtime()-beg_time<<" seconds\n";
		report<<=1;
	}
}

template<typename vertex_t, typename index_t>
bool deep_than_two(
	vertex_t src,
	vertex_t *currq,
	vertex_t *nextq,
	depth_t *depth,
	index_t *adj_card,
	index_t *beg_pos,
	vertex_t *adj_list
){
	
	index_t curr_sz=1;
	index_t next_sz=0;
	currq[0]=src;
	depth[src]=0;
	//sequential check if traversal larger than 2
	for(index_t level=0;;level++)
	{
		next_sz=0;
		for(index_t i=0;i<curr_sz;i++)
		{
			vertex_t frontier=currq[i];
			index_t adj_off=beg_pos[frontier];
			for(index_t j=0;j<adj_card[frontier];j++)
			{
				vertex_t neighbor=adj_list[adj_off+j];
				if(depth[neighbor]==INFTY)
				{
					nextq[next_sz]=neighbor;
					depth[neighbor]=0;
					next_sz++;
				}
			}
		}
		
		if(level==1)
			if(next_sz!=0) return false;
			else return true;
		
		memcpy(currq,nextq,sizeof(vertex_t)*next_sz);
		curr_sz=next_sz;
	}
}


template<	typename vertex_t,
			typename index_t>
void arrange_src
(
	int world_sz,
	int my_id,
	vertex_t *glb_src_list, 
	index_t glb_src_count, 
	vertex_t* &src_list,
	index_t &src_count){
	
	src_count=glb_src_count/world_sz;
	if(glb_src_count != src_count * world_sz){
		std::cout<<"Work not even!\n";
		exit(-1);
	}

	src_list=new vertex_t[src_count];
	int src_off=my_id*src_count;

	memcpy(src_list,glb_src_list+src_off, 
	sizeof(vertex_t)*src_count);
}



template<	typename vertex_t,
			typename index_t>
int read_src_v(	
	vertex_t* 	&glb_src_list,
	index_t			&glb_src_count,
	std::string filename
)
{
	std::stringstream ss;
	std::string str;
	std::ifstream myfile(filename.c_str());
	
	//count num source in input.dat
	//one line is a source
	if(!myfile.good()) return FILE_NOT_EXIST;
	glb_src_count=0;
	while(std::getline(myfile, str))
		glb_src_count++;
	myfile.close();
	//std::cout<<"Read input.dat: "<<glb_src_count<<"\n";

	glb_src_list	= new vertex_t[glb_src_count];
	index_t ptr	= 0;
	myfile.open(filename.c_str());
	while(std::getline(myfile, str)){
		ss.str("");
		ss.clear();
		ss<<str;
		ss>>glb_src_list[ptr];
		ptr++;
	}
	myfile.close();

	return FILE_EXIST;
}

//this function aims to get the vertex list from ss
template <typename vertex_t, typename index_t>
int proc_json(
	std::string kron_addr, 
	index_t			vert_count,
	index_t			edge_count,
	vertex_t*		adj_list,
	index_t*		adj_card,
	index_t*		beg_pos)
{
	
	double beg_time;
	std::stringstream ss;
	std::string		sline;
	index_t i;
	index_t *card_sum		= new index_t[OMP_THDS];
	beg_time=wtime();
	//read card and formulate beg_pos 
	std::cout<<"\n\nReading cardinality ... \n";
	#pragma omp parallel \
	private(i) num_threads(OMP_THDS)
	{
		index_t report=1;
		index_t prc_line=0;
		vertex_t beg_vert,end_vert,vert_id;
		std::string line;
		std::stringstream ss;
		int tid=omp_get_thread_num();
		ss.str("");
		ss.clear();

		//+---------
		//Assuming big file is partioned into small 
		//	files continuously and named sequentially.
		//+----
		//example 1,2,3,4,5 partioned into three files
		//file-0: 1,2; file-1: 3,4; file-2: 5;
		//+---------------------------------------
		ss<<kron_addr<<"/kron."<<tid<<".dat";
		std::string filename=ss.str();
		card_sum[tid]=0;
		std::ifstream file(filename.c_str());
		i=0;beg_vert=0;end_vert=0;
		if(file.is_open()){
			while(std::getline(file, line)){
				progress<index_t>(report,beg_time,tid,prc_line);
				ss.str("");
				ss.clear();
				ss<<line;
				ss>>vert_id;
				ss>>adj_card[vert_id];
				card_sum[tid]+=adj_card[vert_id];
				if(prc_line==0) beg_vert=vert_id;
				prc_line++;
			}
			end_vert=vert_id;
			file.close();
		}else{std::cout<<"card file wrong\n";exit(-4);}
		#pragma omp barrier
		
		beg_pos[beg_vert]=0;
		for(i=0;i<tid;i++)
			beg_pos[beg_vert]	+= card_sum[i];

		for(i=beg_vert+1;i<=end_vert;i++)
			beg_pos[i]=beg_pos[i-1]+adj_card[i-1];
	}
	
	delete[] card_sum;
	std::cout<<"\n\nReading adj_list ... \n";
	#pragma omp parallel \
	private(i) num_threads(OMP_THDS)
	{
		index_t report=1;
		index_t prc_line=0;
		index_t temp;
		vertex_t vert_id;
		std::string line;
		std::stringstream ss;
		int tid=omp_get_thread_num();
		ss.str("");
		ss.clear();
		
		ss<<kron_addr<<"/kron."<<tid<<".dat";
		std::string filename=ss.str();
		std::ifstream file(filename.c_str());
		i=0;
		if(file.is_open()){
			while(std::getline(file, line)){
				progress<index_t>(report,beg_time,tid,prc_line);
				//top-down
				ss.str("");
				ss.clear();
				ss<<line;
				ss>>vert_id;
				ss>>temp;//card
				if(temp!=beg_pos[vert_id+1]-beg_pos[vert_id]){
					std::cout<<tid<<" Wrong 1\n"
					<<"Vert: "<<vert_id<<"\n"
					<<"Expect: "<<beg_pos[vert_id+1]-beg_pos[vert_id]
					<<"\tGot: "<<temp<<"\n";exit(-1);}
				for(i=0;i<beg_pos[vert_id+1]-beg_pos[vert_id];i++)
					ss>>adj_list[beg_pos[vert_id]+i];
				prc_line++;
			}
			file.close();
		}else{std::cout<<"kron file wrong\n";exit(-4);}
	}
	return 0;
}

template <typename vertex_t, typename index_t, typename depth_t>
graph<vertex_t, index_t, depth_t>
::graph(
	graph_reader *g,
		index_t			sw_level,
	index_t			vert_count,
	index_t			edge_count,
	index_t			num_groups,
	index_t			num_agg_bfs,
	int world_sz, int my_id,
	index_t			gpu_id,
	index_t			sml_shed,
	index_t			lrg_shed
):gpu_id(gpu_id),vert_count(vert_count),edge_count(edge_count),
	num_agg_bfs(num_agg_bfs),sml_shed(sml_shed),lrg_shed(lrg_shed),
	sw_level(sw_level), world_sz(world_sz), my_id(my_id)
{
	
	index_t glb_src_count;
	vertex_t *glb_src_list;

	//traversal data
	gdata = new tdata[num_agg_bfs];
	for(index_t i=0;i<num_agg_bfs;i++)
		gdata[i]=new Traverse;
	
	//gpu graph
	ggraph=new CSR_Graph;
	ggraph->vert_count=g->vert_count;
	ggraph->edge_count=g->edge_count;
	ggraph->adj_list=g->csr;
	ggraph->adj_card=new index_t[vert_count];
	ggraph->beg_pos=g->beg_pos;
	EDGES_C = edge_count;//bad programming style, related to direction-switch

	
	for(vertex_t i = 0; i < g->vert_count; i++)
		ggraph->adj_card[i] = g->beg_pos[i+1] - g->beg_pos[i];

	double timer=wtime();
	std::cout<<"Loading edges ... \n";
	//proc_json<vertex_t, index_t>(
	//		kron_addr,
	//		vert_count,
	//		edge_count,
	//		ggraph->adj_list,
	//		ggraph->adj_card,
	//		ggraph->beg_pos);
	std::cout<<"Loading edges time: "<<wtime()-timer<<" second(s)\n";
	std::cout<<"Assuming source vertices are stored in input.dat\n";
	if(true){
		std::cout<<"input.dat is not detected.\n"
		<<"Generating source vertices with rand()%vert_count\n";
	
		vertex_t *currq=new vertex_t[vert_count];
		vertex_t *nextq=new vertex_t[vert_count];
		depth_t *depth=new depth_t[vert_count];
		
		//generate glb_src_list
		glb_src_count = num_groups*num_agg_bfs;
		std::cout<<glb_src_count<<"\n";
		glb_src_list 	= new vertex_t[glb_src_count];
		//std::cout<<"haha\n";
		for(index_t i = 0; i<glb_src_count; i++){
			glb_src_list[i] = rand()%vert_count;
			
			//Not orphan
			//And level more than two
			memset(depth,INFTY,sizeof(depth_t)*vert_count);
			if(deep_than_two<vertex_t,index_t>
				(glb_src_list[i],currq,nextq,depth,
				ggraph->adj_card,ggraph->beg_pos,ggraph->adj_list)){
				--i;continue;
			}

			//Not repeated
			for(index_t j = 0; j< i; j++)
				if(glb_src_list[j] == glb_src_list[i]){
					--i;continue;
				}
		}

		delete[] currq;
		delete[] nextq;
		delete[] depth;
	};
	
	#ifdef GROUPBY
	std::cout<<"before\n";
	vertex_t *glb_new_list=new vertex_t[glb_src_count];
	std::cout<<"after\n";
	GroupBy<vertex_t,index_t>(glb_src_count,glb_src_list,glb_new_list,ggraph);	
	memcpy(glb_src_list,glb_new_list,sizeof(vertex_t)*glb_src_count);
	#endif
	
	arrange_src<vertex_t, index_t>(
	world_sz, my_id, glb_src_list, glb_src_count, src_list, src_count);
	std::cout<<"Global-source-count vs src_count: "<<glb_src_count<<" vs "<<src_count<<"\n";

	//alloc depth and parent array
	if(cudaMallocHost((void **)&depth, sizeof(depth_t)*vert_count))
		std::cout<<"Alloc depth wrong\n";
	parent = new vertex_t[vert_count];

	#ifdef DEBUG_CONSTRUCTOR
	std::cout<<"\n";
	for(int j=0;j<num_lines;j++){
		std::cout<<"Vertex "<<vertex_list[j]->src_ver
			<<"("<<vertex_list[j]->num_conn_ver
			<<","<<vertex_list[j]->depth<<"): ";
		for(int i=0;i<vertex_list[j]->num_conn_ver;i++)
			std::cout<<vertex_list[j]->conn_ver_list[i]<<",";
		std::cout<<"\n";
	}
	#endif

}
