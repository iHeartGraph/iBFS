#ifndef	COMM_HEADER
#define	COMM_HEADER
#include <stdio.h>

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

typedef struct CSR_Graph{
	vertex_t*	adj_list;
	index_t*	adj_card;
	index_t*	beg_pos;
	
	vertex_t*	adj_list_d;
	index_t*	adj_card_d;
	index_t*	beg_pos_d;
	vertex_t*	src_list;	
	
	depth_t*	depth_merge;
	
	//for direction-switching
	vertex_t*	ex_sml_q;
	vertex_t*	ex_mid_q;
	vertex_t*	ex_lrg_q;
	index_t*	ex_cat_sml_sz;
	index_t*	ex_cat_mid_sz;
	index_t*	ex_cat_lrg_sz;
	index_t*	ex_sml_off;
	index_t*	ex_mid_off;
	index_t*	ex_lrg_off;
	
	comp_t		*depth_comp_last;
	comp_t		*depth_comp_curr;
	index_t		vert_count;
	index_t		edge_count;

	cudaStream_t gstream;
} *csr_graph;

typedef struct Traverse{
	depth_t* 	depth_d;
	vertex_t*	parent_d;
	index_t		ex_q_sz;
//	vertex_t*	ex_sml_q;
//	vertex_t*	ex_mid_q;
//	vertex_t*	ex_lrg_q;
	
	index_t	*ex_sml_sz;
	index_t	*ex_mid_sz;
	index_t	*ex_lrg_sz;
	index_t	*ex_sml_sz_d;
	index_t	*ex_mid_sz_d;
	index_t	*ex_lrg_sz_d;
	
///	index_t*	ex_cat_sml_sz;
///	index_t*	ex_cat_mid_sz;
///	index_t*	ex_cat_lrg_sz;
///	index_t*	ex_sml_off;
///	index_t*	ex_mid_off;
///	index_t*	ex_lrg_off;

	index_t*	tr_edges_c_d;
	index_t*	tr_edges_c_h;

	cudaStream_t	*stream;
} *tdata;


#define	VIS				0x02
#define UNVIS			0x00
#define FRT				0x01
#define	SET_VIS(a)		((a)=0x02)
#define	SET_FRT(a)		((a)=0x01)
#define	IS_FRT(a)		((a)==0x01)
#define	IS_VIS(a)		((a)==0x02)
#define	IS_UNVIS(a)		((a)==0x00)

//----------------------------------
//GLOBAL VARIABLES
//---------------------------------
//--------------------------------
#define INFTY			255	
#endif

#ifndef EXTERN
#define EXTERN
//-------------------

index_t ENABLE_BTUP;
index_t SW_LEVEL;
index_t agg_tr_edges;//already traversed edges
index_t EDGES_C;//total edges in the graph

//texture <vertex_t,	1, cudaReadModeElementType> tex_adj_list[4];
//size_t  tex_adj_off[4];

//texture <index_t, 	1, cudaReadModeElementType> tex_card;
//texture <index_t, 	1, cudaReadModeElementType> tex_strt;
//texture <depth_t, 	1, cudaReadModeElementType> tex_depth;
//
//texture <vertex_t,	1, cudaReadModeElementType> tex_sml_exq;
//texture <vertex_t,	1, cudaReadModeElementType> tex_mid_exq;
//texture <vertex_t,	1, cudaReadModeElementType> tex_lrg_exq;

__device__	index_t		hub_vert[16];
__device__ 	depth_t		hub_stat[16];
__device__	index_t		hub_card[16];

#define HUB_SZ			16	
#define HUB_BU_SZ		16//should be 1.25 of HUB_SZ
								//since there is no status
								//array in __shared__ mem
#define HUB_CRITERIA	0	
#define	Q_CARD	3

#endif
