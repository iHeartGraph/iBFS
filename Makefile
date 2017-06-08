EXE=gpu-ibfs

COMMFLAGS=-O3 --compiler-options -Wall -Xptxas -v
CUCC= "$(shell which nvcc)"

CUFLAGS= -arch=sm_35  ${COMMFLAGS}#-Xptxas -dlcm=cg#disable l1 cache
CUFLAGS+= -ccbin=g++ -Xcompiler -fopenmp


MPC    = "$(shell which mpicxx)"
MPCFLAGS  = -Wall -I"$(shell dirname $(CUCC))/../include" -L"$(shell dirname $(CUCC))/../lib64" -lcudart -fopenmp


ifeq ($(enable_monitor), 1)
	CUFLAGS+= -DENABLE_MONITORING
endif

ifeq ($(enable_check), 1)
	CUFLAGS+= -DENABLE_CHECKING
endif


ifeq ($(enable_groupby), 1)
	CUFLAGS+= -DGROUPBY
endif


OBJS=  	main.o \
		ibfs.o \
		reporter.o 

DEPS= 	Makefile \
		expander.cuh \
		inspector.cuh \
		comm.h \
		graph.cuh \
		bfs_gpu_opt.cuh \
		wtime.h \
		validate.h \
		scan.cuh \
		allocator.cuh 

%.o:%.cpp $(DEPS)
	${MPC} -c  ${MPCFLAGS} $< -o $@

%.o:%.cu $(DEPS)
	${CUCC} -c  ${CUFLAGS} $< -o $@

${EXE}:${OBJS}
	${MPC} ${OBJS} $(MPCFLAGS) -o ${EXE}

clean:
	rm -rf *.o ${EXE}
