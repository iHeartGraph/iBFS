#include "scan.cuh"
#include "expander.cuh"
#include "inspector.cuh"
#include "wtime.h"
#include "validate.h"
#include "ibfs.h"

template <typename vertex_t, typename index_t, typename depth_t> 
bool cont_bfs
(
	index_t level, index_t sw_level,
	index_t *qsz,
	index_t *last_ct,
	index_t num_agg_bfs,
	tdata *gdata
){
	bool cont_traverse=false;
	qsz[0]=gdata[0]->ex_sml_sz[0]
					+gdata[0]->ex_mid_sz[0]
					+gdata[0]->ex_lrg_sz[0];
	
	if(!ENABLE_BTUP){
		if(qsz[0]!=0){cont_traverse=true;}
	}else{
		if(last_ct[0]!=qsz[0]){
			cont_traverse=true;
			last_ct[0]=qsz[0];
		}else if((qsz[0]!=0)&&((level==sw_level+1)||(level==sw_level+2)||(level==sw_level+3)||(level==sw_level+4))){
			cont_traverse=true;
			last_ct[0]=qsz[0];
		}
	}
//	std::cout<<"----\n";
//	std::cout<<last_ct[0]<<"\t"<<qsz[0]<<"\n";
//	std::cout<<sw_level<<"\t"<<level<<"\n";
	return cont_traverse;
}

template <typename vertex_t, typename index_t, typename depth_t> 
void ibfs
(  
	vertex_t 	*src_v,
	index_t		src_grp_off,
	tdata			*gdata,
	csr_graph ggraph,
	index_t		num_agg_bfs,
	index_t		vert_count,
	index_t		sw_level,
	index_t		*last_ct,
	depth_t		&level,
	const index_t sml_shed,
	const index_t lrg_shed
)
{
	init_bfs<vertex_t,index_t,depth_t>
	(
		src_v,
		src_grp_off,
		gdata,
		ggraph,
		num_agg_bfs,
		sml_shed,
		lrg_shed
	);
	cudaDeviceSynchronize();

	#ifdef ENABLE_MONITORING
	double tm_insp;
	double tm_expd;
	double tm_step;
	double tm_expand	= 0.0;
	double tm_inspect	= 0.0;
//	depth_t *d_depth;
//	cudaMallocHost((void**)&d_depth,sizeof(depth_t)*vert_count);
//	cudaMemcpy(d_depth,ggraph->depth_merge, sizeof(depth_t)*vert_count,cudaMemcpyDeviceToHost);
//	for(int i=0;i<vert_count;i++)
//		if(d_depth[i]==0)std::cout<<"level 0: "<<i<<"\n";
//	level=0;
//	std::stringstream ss;
//	std::string filename[num_agg_bfs];
//	std::ofstream ofile;
//	for(index_t i=0;i<num_agg_bfs;i++){
//		ss.str("");ss.clear();
//		ss<<"console-"<<src_v[src_grp_off*num_agg_bfs+i]<<".log";filename[i]=ss.str();
//	}
//	for(index_t i=0;i<num_agg_bfs;i++){
//		ofile.open(filename[i].c_str(), std::fstream::out|std::fstream::app);
//		ofile		<<"@Level-"<<(int)level<<"-nextq-insp-expd-step: "
//						<<gdata[0]->ex_sml_sz[0]
//						+gdata[0]->ex_mid_sz[0]
//						+gdata[0]->ex_lrg_sz[0]<<" "
//						<<tm_insp<<" "<<tm_expd<<" "<<tm_step<<"\n";
//		ofile.close();
//	}
//	ofile.close();
	#endif
	
	index_t qsz[num_agg_bfs];
	for(level=1;;level++){
	//	if(level>sw_level+1)
			if(!cont_bfs<vertex_t,index_t,depth_t>
				(level, sw_level, qsz,last_ct,num_agg_bfs,gdata))break;
		
		//----------------------	
		//Expand ex_qs and mark frontier in depth_d	
		//------------------------------------------------
		#ifdef ENABLE_MONITORING
		tm_step=wtime();	
		tm_expd=wtime();	
		#endif
		expander
		<vertex_t, index_t, depth_t,comp_t>
		(
			gdata,
			ggraph,
			num_agg_bfs,
			sw_level,
			level
		);
		cudaDeviceSynchronize();

		#ifdef ENABLE_MONITORING
		tm_expd=wtime()-tm_expd;	
		tm_insp=wtime();	
		#endif
		//-----------------------
		//Generate ex_qs from depth_d
		//---------------------------------
		inspector
		<vertex_t, index_t, depth_t,comp_t>
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
		cudaDeviceSynchronize();
		#ifdef ENABLE_MONITORING
		tm_insp=wtime()-tm_insp;
		tm_step=wtime()-tm_step;
	
    std::cout<<"@Level-"<<(int)level;
    if(level<sw_level) std::cout<<"-td-";
    else if(level==sw_level) std::cout<<"-sw-";
    else std::cout<<"-bu-";
    std::cout<<"nextq-insp-expd-step: "
              <<qsz[0]<<" "<<tm_insp<<" "<<tm_expd<<" "<<tm_step<<"\n";	
	//	for(index_t i=0;i<num_agg_bfs;i++){
	//	//	if(level<=sw_level+1) qsz[i]=0;
	//		ofile.open(filename[i].c_str(), std::fstream::out|std::fstream::app);
	//		ofile		<<"@Level-"<<(int)level<<"-nextq-insp-expd-step: "
	//						<<gdata[0]->ex_sml_sz[0]
	//						+gdata[0]->ex_mid_sz[0]
	//						+gdata[0]->ex_lrg_sz[0]<<" "
	//						<<tm_insp<<" "<<tm_expd<<" "<<tm_step<<"\n";
	//		ofile.close();
	//	}
	//	ofile.close();
		tm_expand		+= tm_expd;
		tm_inspect	+= tm_insp;
		#endif
	}
	#ifdef ENABLE_MONITORING
//	for(index_t i=0;i<num_agg_bfs;i++){
//		ofile.open(filename[i].c_str(), std::fstream::out|std::fstream::app);
//		ofile<<"Total-insp-expd "<<tm_inspect<<" "<<tm_expand<<"\n";
//		ofile.close();
//	}
//	ofile.close();
	std::cout<<"Break-fq-sz: "<<qsz[0]<<"\n";
	std::cout<<"Total-insp-expd "<<tm_inspect<<" "<<tm_expand<<"\n";
	#endif
}

/////////////////////
//CALLING FUNCTION FROM CPU
///////////////////////////
template<typename vertex_t, typename index_t, typename depth_t>
int graph<vertex_t, index_t, depth_t>::
bfs_gpu_coalescing_mem()
{
	cudaSetDevice(gpu_id);
	depth_t *temp;
//	index_t agg_tr_edges, agg_tr_v;
	double tm_bfs;
//	double average_teps	= 0.0;
//	double curr_teps	= 0.0;
//	index_t validate_count = 0;
	
	cudaMallocHost((void **)&temp,sizeof(depth_t)*vert_count*num_agg_bfs);
	for(index_t i=0;i<vert_count*num_agg_bfs;i++) temp[i]=INFTY;
	
	index_t *last_ct=new index_t[num_agg_bfs];
	index_t i;
	for(i=0; i<src_count/num_agg_bfs; i++){
	//	std::cout<<"Test "<<i+1<<"\n";
	//	std::cout<<"Started from: "<<src_list[i]<<"\n";
		ENABLE_BTUP		=false;
		SW_LEVEL			=false;
		agg_tr_edges	=0;
		cudaMemcpy(ggraph->depth_merge,temp,sizeof(depth_t)*vert_count*num_agg_bfs,
					cudaMemcpyHostToDevice);
		cudaMemset(ggraph->depth_comp_last,0,sizeof(comp_t)*vert_count);
		cudaMemset(ggraph->depth_comp_curr,0,sizeof(comp_t)*vert_count);
		last_ct[0]=-1;

		depth_t level=0;
		tm_bfs=wtime();
		ibfs<vertex_t, index_t, depth_t>
		(
			src_list,
			i,
			gdata,
			ggraph,
			num_agg_bfs,
			vert_count,
			sw_level,
			last_ct,
			level,
			sml_shed,
			lrg_shed
		 );
		tm_bfs=wtime()-tm_bfs;
		reporter(tm_bfs, my_id, i);
	//	std::cout<<"Traversal-time: "<<tm_bfs<<" second(s)\n";
		
	//		if(level>2){
	//			validate_count ++;
	//			if(cudaMemcpy(depth, ggraph->depth_merge, 
	//								sizeof(depth_t)*vert_count, 
	//								cudaMemcpyDeviceToHost))
	//				std::cout<<"copy result error\n";

	//			int ret = validate<index_t, vertex_t, depth_t>
	//				(depth, adj_list, adj_card, beg_pos, vert_count);

	//			std::cout<<"\nBFS result validation: "<<
	//						((ret == 0 )? "CORRECT (":"WRONG (")<<ret<<")\n";
	//			report<vertex_t, index_t, depth_t>
	//				(agg_tr_edges, agg_tr_v, depth, adj_card, vert_count);
	//			curr_teps	= agg_tr_edges/(1000000000*tm_bfs);
	//			average_teps= (curr_teps+average_teps*(validate_count-1))
	//								/validate_count;

	//			std::cout<<"Traversed vertices: "<< agg_tr_v<<"\t\t\t"
	//					 <<"Traversed edges: "<<agg_tr_edges<<"\n"
	//					 <<"Traversed time(s) :"<<tm_bfs<<"\t\t"
	//					 <<"Current TEPS (Billion): "<<curr_teps<<"\n"
	//					 <<"Average TEPS (Billion): "<<average_teps<<"\n";
	//		}else{
	//			std::cout<<"Traverse depth is "<<(int)level;
	//		}
	//		std::cout<<"\n====================================\n";
		}
	
//	std::cout<<"Final Average TEPS (Billion): "<<average_teps<<"\n";
	return 0;
}
