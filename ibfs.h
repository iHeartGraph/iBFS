#ifndef IBFSH
#define IBFSH

typedef int index_t;
int ibfs(int world_sz, int my_id, int gpu_count, int args, char *argv[]); 
int reporter(double tm, int my_id, index_t iter);

#endif
