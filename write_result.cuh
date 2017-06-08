#include "graph.h"
#include <stdlib.h>
#include <time.h>
#include <sstream>
template<typename index_t, typename vertex_t, typename depth_t>
int graph<index_t, vertex_t, depth_t>::
write_result()
{
	std::stringstream ss;
	srand(time(NULL));
	std::ofstream result_file;
	std::string file_str="bfs_result.";
	ss<<rand()%8959;
	file_str.append(ss.str());
	file_str.append(".log");
	result_file.open(file_str.c_str());
	
	for(depth_t i=0;;i++){
		index_t counter=0;
		result_file<<"Level "<<(int)i<<"(remaining  vertices) :";
		for(index_t j=0;j<vert_count;j++)
			if(depth[j]==i){
				counter==0 ? :result_file<<",";
				result_file<<j;
				counter++;
			}
		
		result_file<<"----------Total: "<<counter;
		if(!counter) 
			break;
		result_file<<"\n";
	}

	result_file.close();
	return 0;
}
