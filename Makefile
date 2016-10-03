exe=gpu_ibfs.bin

commflags=-O3 --compiler-options -Wall -Xptxas -v
cucc= "$(shell which nvcc)"

cuflags= -arch=sm_35  ${commflags}#-xptxas -dlcm=cg#disable l1 cache
cuflags+= -ccbin=g++ -Xcompiler -fopenmp

ifeq ($(debug), 1)
	cuflags+= -DENABLE_MONITORING
endif

ifeq ($(enable_check), 1)
	cuflags+= -DENABLE_CHECKING
endif


objs = $(patsubst %.cu,%.o,$(wildcard ./*.cu)) \
				$(patsubst %.cpp,%.o,$(wildcard ./*.cpp))

deps = $(wildcard ./*.h) \
			$(wildcard ./*.cuh) \
				Makefile

%.o:%.cu $(deps)
	${cucc} -c ${cuflags} $< -o $@

%.o:%.cpp $(deps)
	${cucc} -c ${cuflags} $< -o $@

${exe}:${objs}
	${cucc} ${objs} $(cuflags) -o ${exe}

clean:
	rm -rf *.o ${exe}
