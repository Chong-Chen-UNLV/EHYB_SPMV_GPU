
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


static void sortRordrList(unsigned int dimension,
		unsigned int nparts,
		unsigned int* part_boundary, 
		unsigned int* rodr_list, 
		unsigned int* numInRow){

	rowS* rodrSVec = (rowS*) malloc(dimension*sizeof(rowS));
	for(unsigned int i = 0; i < dimension; ++i){
		//the location of elements at SVec will be changed 
		//record the original idx since it will be used for
		//new rodr_list
		rodrSVec[rodr_list[i]].idx = i; 
		rodrSVec[rodr_list[i]].nonzeros = numInRow[i];	
	}	
	for(unsigned int i = 0; i < nparts; ++i){
		qsort(&rodrSVec[part_boundary[i]], (part_boundary[i+1] - part_boundary[i]), sizeof(rowS), rowSCompare);	
	}	
	for(unsigned int i = 0; i < dimension; ++i){
		rodr_list[rodrSVec[i].idx] = i;
	}
	free(rodrSVec);
}

/*reorder function with I_rodr, J_rodr, v_rodr, rodr_list as output*/
void matrixReorder(matrixCOO* localMatrixCOO,
		int* reorderList,
		cb_s cb)
{
	int dimension = inputMatrix->dimension;
	inputMatrix->partBoundary = (int*) malloc(dimension*sizeof(int));
	localMatrixCOO->nparts = ceil(dimension/vectorCacheSize);
	int* partBoundary = inputMatrix->partBoundary;
	int nparts = localMatrixCOO->nparts;
	printf("nparts is %d\n", nParts);
	int tempI, tempIdx, tempJ;
	int maxCol;
	double tempV;
	int* I = localMatrixCOO->I;
	int* J = localMatrixCOO-J;
	int* V= localMatrixCOO->V;
	int* rowIdx = localMatrixCOO->rowIdx; 
	int* reorderList = localMatrixCOO->reorderList;
	
	int* newI = (int*)malloc(sizeof(int)*localMatrixCOO->totalNum);
	int* newJ = (int*)malloc(sizeof(int)*localMatrixCOO->totalNum);
	double* newV = (double*)malloc(sizeof(double)*localMatrixCOO->totalNum);

	int* numInRow2 = localMatrixCOO->numInRow2;
	/*transfer the COO format to CSR format, do the partitioning*/
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
			&nparts,
			NULL,
			&ubvec,
			options,
			&r_edgecut,
			partVec);
	/*do the reordering based on partition*/
	printf("partition finished\n");
	gettimeofday(&end, NULL);
	
	printf("partition time is %ld us\n",(end.tv_sec * 1000000 + end.tv_usec)-(start.tv_sec * 1000000 + start.tv_usec));
	int numInRow2[dimension] = 0;	
	int* partFilled = (int* )calloc(nparts, sizeof(int)); 	
	for(int i = 0; i < dimension; ++i){
		partSize[partVec[i]] += 1;
	}
	partBias[0] = 0;
	for(int i = 1; i < nparts + 1; ++i){
		partBias[i] = partBias[i-1] + partSize[i-1];
	}
	int permIdx;	
	for(int i = 0; i < dimension; i++){
		permIdx = partFilled[partVec[i]] + partBias[partVec[i]];
		reorderList[i] = permIdx;
		partFilled[partVec[i]]+=1;
		numInRow2[i] = 0;
	}	
	 
	for(int i = 0; i <= nparts; i++){
		inputMatrix->partBoundary[i] = partBias[i];
	}
	//numInRow2 means "real" numInRow in the  
	//block ELL part of EHYB, notice that some of the
	//elements will still go to ER part even its col 
	//index belongs to same block of this row 
	for(int i = 0; i < totalNum; i++){
		if(partVec[J[i]] == partVec[I[i]]){
			numInRow2[I[i]]+=1;
		}
	}
	//sort the reorder list by number of nozero per row
	//it may not working well on factorized preconditioner
	if(cb.SORT){
		sortRordrList(dimension, nparts, partBoundary, reorderList, numInRow2);
	}
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
	for(int i = 0; i < localMatrixCOO->totalNum; i++){
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
	localMatrixCOO->I=newI;
	localMatrixCOO->J=newJ;
	localMatrixCOO->V=newV;
	free(cutSize);
	free(partSize);
	free(partVec);
	free(cwghts);
	free(partBias);
	free(partFilled);
}

void vectorReorder(const unsigned int dimension, const double* v_in, 
			double* v_rodr, const unsigned int* rodr_list){
	for(int i=0; i < dimension; i++) v_rodr[rodr_list[i]] = v_in[i];	
}
void vector_recover(const unsigned int dimension, const double* v_rodr, double* v, const unsigned int* rodr_list){
	unsigned int* rodr_list_recover= (unsigned int*) malloc(dimension*sizeof(unsigned int));
	for(int i=0; i < dimension; i++) rodr_list_recover[rodr_list[i]] = i;	
	for(int i=0; i < dimension; i++) v[rodr_list_recover[i]] = v_rodr[i];	
	free(rodr_list_recover);
}

void update_numInRowL(const unsigned int totalNum, 
			const unsigned int dimension, 
			unsigned int* I_rodr, 
			unsigned int* J_rodr, 
			double* V_rodr, 
			unsigned int* numInRowL,
			unsigned int* maxRowNumL,
			unsigned int* maxRowNumLP,
			unsigned int* row_idxL, 
			unsigned int* row_idxLP,
			double* diag){
	
	unsigned int* numInRowLP_ = (unsigned int*)malloc(dimension*sizeof(unsigned int));
	unsigned int maxL = 0;
	unsigned int maxLP = 0;
	for(unsigned int i = 0; i < dimension; ++i){
		numInRowL[i] = 0;		
		numInRowLP_[i] = 0;		
		row_idxL[i] = 0;
		row_idxLP[i] = 0;
	}
	row_idxL[dimension] = 0;
	row_idxLP[dimension] = 0;

	for(unsigned int i=0; i< totalNum; i++){
		if(I_rodr[i] == J_rodr[i]){
			diag[I_rodr[i]] = V_rodr[i];
			numInRowL[I_rodr[i]] += 1;	
			numInRowLP_[I_rodr[i]] += 1;	
		}
		else if(I_rodr[i] < J_rodr[i]){
			numInRowL[I_rodr[i]] += 1;	
		}
		else{
			numInRowLP_[I_rodr[i]] += 1;	
		}
	}
	for(unsigned int i = 1; i <= dimension; ++i){
		row_idxL[i]=row_idxL[i-1]+numInRowL[i-1];
		row_idxLP[i]=row_idxLP[i-1]+numInRowLP_[i-1];
		if(numInRowLP_[i-1] > maxLP)
			maxLP = numInRowLP_[i-1];
		if(numInRowL[i-1] > maxL)
			maxL = numInRowL[i-1];
	}
	*maxRowNumLP = maxLP;
	*maxRowNumL = maxL;
	free(numInRowLP_);
}
