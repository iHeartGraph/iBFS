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

#ifndef __UTIL_H__
#define	__UTIL_H__
#include <stdint.h>
#include <sys/stat.h>
#include <stdio.h>

inline off_t fsize(const char *filename) {
	struct stat st; 
	if (stat(filename, &st) == 0)
		return st.st_size;
	return -1; 
}

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", \
        cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define H_ERR( err ) \
  (HandleError( err, __FILE__, __LINE__ ))

//////////////////////////////////////////////////
//SCALE*THDS_NUMS*sizeof(int) should be 
//limited by the size of the shared memory
/////////////////////////////////////////////////
//////////////////////////////////////////////////
#define	THDS_NUM			256	
#define	BLKS_NUM			256	
#define OMP_THDS			24
#define	V_NON_INC			-1
#define V_INI_HUB			-2

#define VALIDATE_TIMES		1
#define NUM_SRC		 		128	
#define INFTY					255	
#define	Q_CARD				3
enum ex_q_t
{
  SML_Q,
  MID_Q,
  LRG_Q,
  NONE
};

//--------------------------
typedef 	int 						index_t;
typedef		int							vertex_t;
typedef 	unsigned char 	depth_t;
typedef 	uint4	comp_t;
//--------------------------------

//operator overloading
inline __host__ __device__ uint4 operator-(uint4 a, uint4 b)
{
  return make_uint4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}

inline __host__ __device__ uint4 operator|=(uint4 a, uint4 b)
{
  return make_uint4(a.x |= b.x, a.y |= b.y, a.z |= b.z,  a.w |= b.w);
}

inline __host__ __device__ uint4 operator|(uint4 a, uint4 b)
{
  return make_uint4(a.x | b.x, a.y | b.y, a.z | b.z,  a.w | b.w);
}

inline __host__ __device__ uint4 operator^(uint4 a, uint4 b)
{
  return make_uint4(a.x ^ b.x, a.y ^ b.y, a.z ^ b.z,  a.w ^ b.w);
}

inline __host__ __device__ uint4 operator~(uint4 a)
{
  return make_uint4(~a.x, ~a.y, ~a.z, ~a.w);
}

inline __host__ __device__ uint4 make_uint4(int a)
{
  return make_uint4(a, 0, 0, 0);
}

inline __host__ __device__ uint4 operator<<(uint4 a, int b)
{
  if(b<32) a.x=(a.x<<b);
	else
	{
		if(b<64) a.y=a.x<<(b-32);
		else if(b<96) a.z=a.x<<(b-64);
		else a.w=a.x<<(b-96);
		
		a.x=0;
	}

	return a;
}

inline __host__ __device__ uint4 operator&(uint4 a, uint4 b)
{
	return make_uint4(a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w);
}

inline __host__ __device__ int operator!=(uint4 a, unsigned int b)
{
	if(	a.x != b ||
			a.y != b ||
			a.z != b ||
			a.w != b ) return true;
	else return false;
}


inline __host__ __device__ int operator!=(uint4 a, uint4 b)
{
	if(	a.x != b.x ||
			a.y != b.y ||
			a.z != b.z ||
			a.w != b.w ) return true;
	else return false;
}

inline __host__ __device__ int operator==(uint4 a, unsigned int b)
{
	if(a.x == b && 
		 a.y == b &&
		 a.z == b &&
		 a.w == b) return true;
	else return false;
}

#endif
