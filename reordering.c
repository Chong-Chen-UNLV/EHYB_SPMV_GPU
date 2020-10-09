
/*
 * Reordering functions, using paritioning libraries to conduct
 * graphic partition, then the parition result will be used
 * for reordering. This reordering will divide matrices into
 * different parts, and reduce the dependent between parts 
 * when conduct iterative SpMV operations.
 * The library used (or will be used) here includes mu-metis,
 * PaToH (third part library), and in to future I will use 
 * graphic partition methods developed by me 
 * */

#include "Partition.h"
#include "solver.h"
#include "kernel.h"

#include<string.h>


static void sortRordrList(int dimension,
		int nParts,
		int* part_boundary, 
		int* rodr_list, 
		int* numInRow){

	rowS* rodrSVec = (rowS*) malloc(dimension*sizeof(rowS));
	for(unsigned int i = 0; i < dimension; ++i){
		//the location of elements at SVec will be changed 
		//record the original idx since it will be used for
		//new rodr_list
		rodrSVec[rodr_list[i]].idx = i; 
		rodrSVec[rodr_list[i]].nonzeros = numInRow[i];	
	}	
	for(unsigned int i = 0; i < nParts; ++i){
		qsort(&rodrSVec[part_boundary[i]], (part_boundary[i+1] - part_boundary[i]), sizeof(rowS), rowSCompare);	
	}	
	for(unsigned int i = 0; i < dimension; ++i){
		rodr_list[rodrSVec[i].idx] = i;
	}
	free(rodrSVec);
}

/*reorder function with I_rodr, J_rodr, v_rodr, rodr_list as output*/
void matrixReorder(matrixCOO* inputMatrix)
{
	unsigned int dimension = inputMatrix->dimension;
	inputMatrix->partBoundary = (int*) malloc(dimension*sizeof(int));
	inputMatrix->nParts = ceil(dimension/vectorCacheSize);
	int* partBoundary = inputMatrix->partBoundary;
	unsigned int nParts = inputMatrix->nParts;
	printf("nParts is %d\n", nParts);
	int tempI, tempIdx, tempJ;
	int maxCol;
	double tempV;
	uint32_t* J = (uint32_t*)malloc(sizeof(uint32_t)*inputMatrix->totalNum);
	for(int i = 0; i < inputMatrix->totalNum; ++i){
		J[i] = inputMatrix->J[i];
	}
	free(inputMatrix->J);
	int* I = inputMatrix->I;
	double* V= inputMatrix->V;
	uint32_t* rowIdx = (uint32_t*)malloc(sizeof(uint32_t)*dimension); 
	int* reorderList = inputMatrix->reorderList;
	int* numInRow = inputMatrix->numInRow;
	
	int* newI = (int*)malloc(sizeof(uint32_t)*inputMatrix->totalNum);
	int* newJ = (int*)malloc(sizeof(int)*inputMatrix->totalNum);
	double* newV = (double*)malloc(sizeof(double)*inputMatrix->totalNum);

	int* numInRow2 = inputMatrix->numInRow2;
	/*transfer the COO format to CSR format, do the partitioning*/
	unsigned int *partVec, *cwghts;
	partVec = (unsigned int *) calloc(dimension, sizeof(unsigned int));
	cwghts = (unsigned int *) calloc(dimension, sizeof(int));

	for(int i=0; i < dimension; i++){
		cwghts[i] = 1;
		rowIdx[i] = inputMatrix->rowIdx[i]; 
	}
	partVec = (unsigned int *) calloc(dimension, sizeof(int));
	//partweights = (int *)calloc(nParts*sizeof(int));
	int* cutSize = (int *)calloc(nParts, sizeof(int));
	int* partSize = (int *)calloc(nParts, sizeof(int));
	int* partBias = (int *)calloc(nParts + 1, sizeof(int));
	double* options = mtmetis_init_options();
	
	unsigned int ncon = 1;
	float ubvec = 1.001;
	options[MTMETIS_OPTION_NTHREADS] = 16;
	mtmetis_wgt_type r_edgecut;
	struct timeval start, end;
	gettimeofday(&start, NULL);
	/*call the graph patition library*/
	printf("start k-way partition\n");
	MTMETIS_PartGraphKway(
			&dimension,
			&ncon,
			rowIdx,
			J,
			NULL,
			NULL,
			NULL,
			&nParts,
			NULL,
			&ubvec,
			options,
			&r_edgecut,
			partVec);
	/*do the reordering based on partition*/
	printf("partition finished\n");
	gettimeofday(&end, NULL);
	
	printf("partition time is %ld us\n",(end.tv_sec * 1000000 + end.tv_usec)-(start.tv_sec * 1000000 + start.tv_usec));
	numInRow2[dimension] = 0;	
	int* partFilled = (int* )calloc(nParts, sizeof(int)); 	
	for(int i = 0; i < dimension; ++i){
		partSize[partVec[i]] += 1;
	}
	partBias[0] = 0;
	for(int i = 1; i < nParts + 1; ++i){
		partBias[i] = partBias[i-1] + partSize[i-1];
	}
	int permIdx;	
	for(int i = 0; i < dimension; i++){
		permIdx = partFilled[partVec[i]] + partBias[partVec[i]];
		reorderList[i] = permIdx;
		partFilled[partVec[i]]+=1;
		numInRow2[i] = 0;
	}	
	 
	for(int i = 0; i <= nParts; i++){
		inputMatrix->partBoundary[i] = partBias[i];
	}
	//numInRow2 means "real" numInRow in the  
	//block ELL part of EHYB, notice that some of the
	//elements will still go to ER part even its col 
	//index belongs to same block of this row 
	for(int i = 0; i < inputMatrix->totalNum; i++){
		if(partVec[J[i]] == partVec[I[i]]){
			numInRow2[I[i]]+=1;
		}
	}
	sortRordrList(dimension, nParts, partBoundary, reorderList, numInRow2);
	//reordering need two iteration
	for(int i = 0; i < dimension ; i++){
		rowIdx[reorderList[i]] = rowIdx[i-1] + numInRow[i-1];
	}
	memset(numInRow, 0 , dimension*sizeof(int));
	rowIdx[0] = 0;
	for(int i= 1; i<=dimension; i++){
		rowIdx[i] = rowIdx[i-1] + numInRow[i-1];
		numInRow[i-1] = 0;
	}
	for(int i = 0; i < inputMatrix->totalNum; i++){
		tempI = reorderList[I[i]];
		tempJ = reorderList[J[i]];	
		tempIdx = rowIdx[tempI] + numInRow[tempI];
		newI[tempIdx] = tempI;
		newJ[tempIdx] = tempJ;
		if(tempI == tempJ && V[i] == 0)
			printf("error happend tempIdx is %d with tempI %d  V is %f\n", tempIdx, tempI, V[i]);
		newV[tempIdx] = V[i];
		numInRow[tempI] += 1;
	}	
	free(I);
	free(J);
	free(V);
	inputMatrix->I=newI;
	inputMatrix->J=newJ;
	inputMatrix->V=newV;
	free(cutSize);
	free(partSize);
	free(partVec);
	free(cwghts);
	free(partBias);
	free(partFilled);
}

void vectorReorder(int dimension, const double* v_in, 
			double* v_rodr, const int* rodr_list)
{
	for(int i=0; i < dimension; i++) v_rodr[rodr_list[i]] = v_in[i];	
}

void vectorRecover(const int dimension, const double* v_rodr, double* v, const int* rodr_list){
	int* rodr_list_recover= (int*) malloc(dimension*sizeof(int));
	for(int i=0; i < dimension; i++) rodr_list_recover[rodr_list[i]] = i;	
	for(int i=0; i < dimension; i++) v[rodr_list_recover[i]] = v_rodr[i];	
	free(rodr_list_recover);
}


