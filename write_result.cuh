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
#include <stdlib.h>
#include <time.h>
#include <sstream>
int gpu_ibfs::
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
				result_file<<(counter==0 ? "":",");
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
