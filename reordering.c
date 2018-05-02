
/*
 * Reordering functions, using paritioning libraries to conduct
 * graphic partition, then the parition result will be used
 * for reordering. This reordering will divide matrices into
 * different parts, and reduce the dependent between parts 
 * when conduct iterative SpMV operations.
 * The library used (or will be used) here includes mu-metis,
 * PaToH, and other graphic partition methods developed by 
 * myself
 * */

#include "Partition.h"

/*reorder function with I_rodr, J_rodr, v_rodr, rodr_list as output*/
void matrix_reorder(const int* dimension_in, const int totalNum, const int* I, const int* J, const float* V, 
		int* numInRow, int* I_rodr, int* J_rodr, float* V_rodr, int* rodr_list, int* part_boundary,
		const int rodr_option, const unsigned int nparts){

	unsigned int dimension = *dimension_in;
	unsigned int *colVal = (unsigned int *) malloc(totalNum*sizeof(int));
	unsigned int* row_idx= (unsigned int*) malloc((dimension + 1)*sizeof(int));
	int tempI, tempIdx, tempJ;
	float tempV;
	/*transfer the COO format to CSR format, do the partitioning*/
	for(int i= 1; i <= dimension; i++){
		row_idx[i] = row_idx[i-1] + numInRow[i-1];
		numInRow[i-1] = 0;
	}
	for(int i =0; i < dimension; i++){
		int rowOfst = row_idx[I[i]] + numInRow[I[i]];
		colVal[rowOfst] = J[i]; 
		numInRow[I[i]] += 1;
	}
	
	unsigned int *partVec, *cwghts;
	partVec = (unsigned int *) calloc(dimension, sizeof(unsigned int));
	cwghts = (unsigned int *) calloc(dimension, sizeof(int));

	for(int i=0; i < dimension; i++){
		cwghts[i] = 1;
	}
	partVec = (unsigned int *) calloc(dimension, sizeof(int));
	//partweights = (int *)calloc(nparts*sizeof(int));
	int* cutSize = (int *)calloc(nparts, sizeof(int));
	int* partSize = (int *)calloc(nparts, sizeof(int));
	int* partBias = (int *)calloc(nparts + 1, sizeof(int));
	double* options = mtmetis_init_options();
	
	int ncon = 1;
	float ubvec = 1;
	options[MTMETIS_OPTION_NTHREADS] = 16;
	mtmetis_wgt_type r_edgecut;
	struct timeval start, end;
	gettimeofday(&start, NULL);
	/*call the graph patition library*/
	//printf("start k-way partition\n");
	MTMETIS_PartGraphKway(
			&dimension,
			NULL,
			row_idx,
			colVal,
			NULL,
			NULL,
			NULL,
			&nparts,
			NULL,
			&ubvec,
			options,
			&r_edgecut,
			partVec);
	/*do the reordering based on partition*/
	int* part_filled = (int* )calloc(nparts, sizeof(int)); 	
	for(int i = 0; i < dimension; ++i){
		partSize[partVec[i]] += 1;
	}
	for(int i = 1; i < nparts + 1; ++i){
		partBias[i] = partBias[i-1] + partSize[i];
	}
	int perm_idx;	
	for(int i=0; i < dimension; i++){
		perm_idx = part_filled[partVec[i]] + partBias[partVec[i]];
		rodr_list[i] = perm_idx;
		part_filled[partVec[i]]+=1;
	}	
	part_boundary[0] = 0;
	for(int i=1; i <= nparts; i++){
		part_boundary[i] = part_boundary[i-1] + part_filled[i];
	}

	for(int i = 0; i < totalNum; i++){
		tempI = rodr_list[I[i]];
		tempJ = rodr_list[J[i]];	
		tempIdx = row_idx[tempI] + numInRow[tempI];
		I_rodr[tempIdx] = tempI;
		J_rodr[tempIdx] = tempJ;
		V_rodr[tempIdx] = V[i];
		numInRow[tempI] += 1;
	}	
	free(colVal);
	free(row_idx);
	free(cutSize);
	free(partSize);
	free(partVec);
	free(cwghts);
	free(partBias);
	free(part_filled);
}

void vector_reorder(const int dimension, const int* rodr_list, const float* v_in, float* v_rodr){
	for(int i=0; i < dimension; i++) v_rodr[rodr_list[i]] = v_in[i];	
}
void vector_recover(const int dimension, const int* rodr_list, const float* v_rodr, float* v){
	int* rodr_list_recover= (int*) malloc(dimension*sizeof(int));
	for(int i=0; i < dimension; i++) rodr_list_recover[rodr_list[i]] = i;	
	for(int i=0; i < dimension; i++) v[rodr_list_recover[i]] = v[i];	
	free(rodr_list_recover);
}

void update_numInRow(const int totalNum, const int dimension, int* I_rodr, 
		int* J_rodr, float* V_odr, int* numInRow, int* numInRowL,
		int* row_idx, int* row_idxL, float* diag){
	
	for(int i = 0; i < dimension; ++i){
		numInRow[i] = 0;		
		row_idx[i] = 0;
		row_idxL[i] = 0;
	}
	row_idx[dimension] = 0;
	row_idxL[dimension] = 0;

	for(int i=0; i< totalNum; i++){
		numInRow[I_rodr[i]] += 1;	
		if(I_rodr[i] == J_rodr[i]){
			diag[I_rodr[i]] = V_odr[i];
			numInRowL[I_rodr[i]] += 1;	
		}
		else if(I_rodr[i] < J_rodr[i]){
			numInRowL[I_rodr[i]] += 1;	
		}
	}
	for(int i = 0; i < dimension; ++i){
		row_idx[i]=row_idx[i-1]+numInRow[i-1];
		row_idxL[i]=row_idxL[i-1]+numInRowL[i-1];
	}
}
