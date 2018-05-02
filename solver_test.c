#include "fspai.h"
#include "mmio.h"
#include "kernel.h"
#include "reordering.h"
#include "solver.h"
#include <omp.h>

#define GPU 0
#define RODR false
#define MAXthread 2 
#define MAXthread2 6

int *I_precond;
int *J_precond;
float *V_precond;

void fspaiCPU(S *SInput);
void fspai(S *SInput);	



int main(int argc, char* argv[])
{	
	int totalNum=0;
	int lowerNum=0;
	int lowerNumPrecond=0;
	int totalNumPrecond=0;
	int totalNumPrecondP=0;
	int ret_code;
	MM_typecode matcode,matcode2;
	int MAXIter;
	FILE *f,*f2;

	int dimension, N, bandwidth,bandwidthPrecond;   
	int maxRowNum, maxRowNumPrecond, maxRowNumPrecondP;
	int *I, *J;
	int *lowerI, *lowerJ;
	float *V;
	float *lowerV;
	float *x;
	float *y;
	float *x_compare;
	float *error_track;
	float result_error;
	int finish;
	char fileName[100];
	char fileName2[100];

	if(argc!=3) {  /* Check command line inputs */
		printf("Usage: distribute [matrixFile]\n");
		exit(0);
	}		
	sprintf(fileName, "../read/%s.mtx", argv[1]);		
	MAXIter=atoi(argv[2]);
	//---------------------------------read the matrix---------------------------
	if ((f = fopen(fileName, "r")) == NULL) 
		exit(1);

	if (mm_read_banner(f, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
			mm_is_sparse(matcode) )
	{
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	/* find out size of sparse matrix .... */

	if ((ret_code = mm_read_mtx_crd_size(f, &dimension, &N, &lowerNum)) !=0)
		exit(1);	

	/*The overall number of nozeros in this matrix*/
	totalNum=lowerNum*2-dimension;
	totalNumPrecond=lowerNum;
	totalNumPrecondP = lowerNum;

	size_t size0=lowerNum*sizeof(int);
	size_t size1=lowerNum*sizeof(float);
	size_t size2=totalNum*sizeof(int);
	size_t size3=totalNum*sizeof(float);
	size_t size4=dimension*sizeof(int);
	size_t size5=dimension*sizeof(float);
	size_t size6=lowerNum*sizeof(int);
	size_t size7=lowerNum*sizeof(float);


	lowerJ=(int *) malloc(size0);
	lowerI=(int *) malloc(size0);
	lowerV=(float *) malloc(size1);
	I=(int *) malloc(size2);
	J=(int *) malloc(size2);
	V=(float *) malloc(size3);
	x=(float *) malloc(size5);
	y=(float *) malloc(size5);
	float *diag=(float *) malloc(size5);
	x_compare=(float *) malloc(size5);
	error_track=(float *) malloc(MAXIter*sizeof(float));

	int *numInRowL;
	int *row_idxL;
	int *numInRowLP;
	int *row_idxLP;	
	int *numInRow;

	numInRowL=(int *) malloc(size4);
	numInRowLP=(int *) malloc(size4);
	row_idxL=(int *) malloc(size4 + sizeof(int));
	row_idxLP=(int *) malloc(size4 + sizeof(int));
	numInRow=(int *) malloc(size4);
	int tempI, tempJ;
	float tempV;
	for (int i=0; i<lowerNum; i++)
	{
		fscanf(f, "%d %d %f\n", &tempI, &tempJ, &tempV);
		lowerJ[i]=tempJ-1;  /* adjust from 1-based to 0-based */
		lowerI[i]=tempI-1;
		lowerV[i]=tempV;
		numInRow[tempI-1]+=1;
		numInRowL[tempJ-1]+=1;
		numInRowLP[tempI-1]+=1;
		if (tempI!=tempJ)
		{
			numInRow[tempJ-1]+=1;
		}		
		if (lowerI[i]-lowerJ[i]>bandwidth) bandwidth=lowerI[i]-lowerJ[i];
	}


	int *row_idx=(int *)malloc((dimension+1)*sizeof(int));
	maxRowNum=0;
	maxRowNumPrecond=0;
	maxRowNumPrecondP=0;
	row_idx[0] = 0;
	row_idxL[0] = 0;
	row_idxLP[0] = 0;
	for (int i=1;i<= dimension;i++)
	{

		if (numInRow[i-1]>maxRowNum)
			maxRowNum=numInRow[i-1];
		if (numInRowL[i-1]>maxRowNumPrecond)
			maxRowNumPrecond=numInRowL[i-1];
		if (numInRowLP[i-1]>maxRowNumPrecondP)
			maxRowNumPrecondP=numInRowLP[i-1];			

		row_idx[i]=row_idx[i-1]+numInRow[i-1];
		row_idxL[i]=row_idxL[i-1]+numInRowL[i-1];
		row_idxLP[i]=row_idxLP[i-1]+numInRowLP[i-1];
		numInRow[i-1]=0;
		numInRowLP[i-1]=0;
		//determine y
	}	
	if (numInRow[dimension-1]>maxRowNum) maxRowNum=numInRow[dimension-1];
	if (numInRowL[dimension-1]>maxRowNumPrecond) maxRowNumPrecond=numInRowL[dimension-1];	
	if (numInRowLP[dimension-1]>maxRowNumPrecondP) maxRowNumPrecondP=numInRowLP[dimension-1];
	numInRow[dimension-1]=0;
	numInRowLP[dimension-1]=0;
	MAXIter=atoi(argv[2]);
	for (int i=0;i<dimension;i++)
	{		
		srand(i);
		x_compare[i]=(float) (rand()%200-100)/100;
		//x_compare[i]=1;
	}
	int index1, index2;

	for (int i=0;i<lowerNum;i++)
	{
		tempI=lowerI[i];
		tempJ=lowerJ[i];
		tempV=lowerV[i];
		index1=row_idx[tempI]+numInRow[tempI];
		index2=row_idx[tempJ]+numInRow[tempJ];
		numInRow[tempI]+=1;
		I[index1]=tempI;
		J[index1]=tempJ;
		V[index1]=tempV;
		y[tempI]+=tempV*x_compare[tempJ];
		if (tempI != tempJ)
		{
			numInRow[tempJ]+=1;
			I[index2]=tempJ;
			J[index2]=tempI;
			V[index2]=tempV;
			y[tempJ]+=tempV*x_compare[tempI];
		}
		else
		{
			diag[tempI]=tempV;
		}
	}	
	/*-----------------do the reordering with metis/hmetis, determine the value------------*/
	/*suffix _rodr means reordered*/		
	int* I_rodr, J_rodr, part_boundary, rodr_list;
	float* V_rodr, x_rodr, y_rodr;
	
	if(RODR){
		rodr_list = (float* )calloc(dimension, sizeof(int)); 
		part_boundary = (int* )calloc((blocks + 1), sizeof(int)); 	
		I_rodr = (int *) malloc(size2);
		J_rodr = (int *) malloc(size2);
		V_rodr = (float *) malloc(size3);
		x_rodr = (float* )calloc(dimension, sizeof(int)); 
		y_rodr = (float* )calloc(dimension, sizeof(int)); 
		matrix_reorder(dimension, totalNum, I, J, V, numInRow, I_rodr, J_rodr, V_rodr, rodr_list, part_boundary, option);
		vector_reorder(dimension, y, y_rodr, rodr_list);
		update_numInRow(totalNum, dimension, I_rodr, J_roder, V_rodr, numInRow, numInRowL, row_idx, row_idxL, diag);
	}
	/*---------------------read the preconditioner ------------------------------*/


	int *I_precond=(int *) malloc(size6);
	int *J_precond=(int *) malloc(size6);
	float *V_precond=(float *) malloc(size7);	

	/*int rt=pthread_barrier_init(&barr, NULL, MAXthread);
	  rt=pthread_barrier_init(&barr1, NULL, MAXthread);*/

	S Sthread[MAXthread];

	for (int t=0;t<MAXthread;t++)
	{
		if(RODR){
			Sthread[t].I=I_rodr;
			Sthread[t].J=J_rodr;
			Sthread[t].V=V_rodr;
		}
		else{
			Sthread[t].I=I;
			Sthread[t].J=J;
			Sthread[t].V=V;
		}
		Sthread[t].I_precond=I_precond;
		Sthread[t].J_precond=J_precond;
		Sthread[t].V_precond=V_precond;
		Sthread[t].maxRowNum=maxRowNum;
		Sthread[t].numInRow=numInRow;
		Sthread[t].row_idx=row_idx;
		Sthread[t].numInRowPrecond=numInRowL;
		Sthread[t].row_idxPrecond=row_idxL;
		Sthread[t].diag=diag;
		if (t==0) Sthread[t].colStart=0;
		else Sthread[t].colStart=(dimension/MAXthread)*t;
		if (t==MAXthread-1) Sthread[t].colEnd=dimension;
		else Sthread[t].colEnd=(dimension/MAXthread)*(t+1);
		Sthread[t].id=t;
	}	
	struct timeval start1, end1;
	struct timeval start2, end2;

	gettimeofday(&start1, NULL);
	omp_set_num_threads(MAXthread);
	int rank;
#pragma omp parallel private(rank)
	{
		rank=omp_get_thread_num();
		fspaiCPU(&Sthread[rank]);
	}
	gettimeofday(&end1, NULL);	
	printf("fspai CPU time is %ld us\n",(end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec));
	int *I_precondP=(int *) malloc(size6);
	int *J_precondP=(int *) malloc(size6);
	float *V_precondP=(float *) malloc(size7);
	
	for (int i=0; i<totalNumPrecond; i++)
	{
		tempI=J_precond[i];
		tempJ=I_precond[i];
		tempV=V_precond[i];
		if (tempI<0||tempI>dimension-1||tempJ<0||tempJ>dimension-1)
			 printf("error happend at %d with tempI %d and tempJ %d\n", i ,tempI, tempJ); 
		index1=row_idxLP[tempI]+numInRowLP[tempI];
		I_precondP[index1]=tempI;
		J_precondP[index1]=tempJ;
		V_precondP[index1]=tempV;
		numInRowLP[tempI]+=1;
	}


	int *I_precond_accum, *I_precondP_accum;
	I_precond_accum = (int *)malloc((dimension + 1)*sizeof(int));
	I_precondP_accum = (int *)malloc((dimension + 1)*sizeof(int));

	int realIter, pp, procNum;
	struct timeval start, end;

#pragma omp parallel
	{
#pragma omp master
		{
			procNum=omp_get_num_threads();
		}
	}
	printf("thread num is %d\n",procNum);
	if(GPU){
		if(RODR){
			solverGPU_HYB(dimension, totalNum, row_idx, J_rodr, V_rodr, totalNumPrecond, 
					I_precond_accum, J_precond, V_precond,totalNumPrecondP, 
					I_precondP_accum, J_precondP, V_precondP, 
					y_rodr, x_rodr, MAXIter, realIter, true, partition_boundary);
		}
		else{
			solverGPU_HYB(dimension, totalNum, row_idx, J, V, totalNumPrecond, 
					I_precond_accum, J_precond, V_precond,totalNumPrecondP, 
					I_precondP_accum, J_precondP, V_precondP, 
					y, x, MAXIter, realIter, false, NULL);
		}
	}
	else{
		if(RODR){
			solverPrecondCPU(dimension, totalNum, row_idx, J_rodr, V_rodr, totalNumPrecond, 
					I_precond_accum, J_precond, V_precond,totalNumPrecondP, 
					I_precondP_accum, J_precondP, V_precondP, 
					y_rodr, x_rodr, MAXIter, realIter);
		}
		else{
			solverPrecondCPU(dimension, totalNum, row_idx, J, V, totalNumPrecond, 
					I_precond_accum, J_precond, V_precond,totalNumPrecondP, 
					I_precondP_accum, J_precondP, V_precondP, 
					y, x, MAXIter, realIter);
		}
	}
	//solverPrecondGPU_CUSPARSE(dimension, totalNum, *row_idx, *J, *V, totalNumPrecond, *I_precond_accum, *J_precond, *V_precond, 
	//	totalNumPrecondP, *I_precondP_accum, *J_precondP, *V_precondP, *y, *x, MAXIter, *realIter);
		
	if(RODR){
		vector_recover(dimension, rodr_list, x_rodr, x);
	}

	for (int i=0;i<10;i++)
	{
		//printf("Xeon_phi I is %d J %d is V is %f\n",I_precond[i+10000], J_precond[i+10000], V_precond[i+10000]);
		//printf("CPU I is %d, J is %d, V is %f\n",I_precond2[i+10000],J_precond2[i+10000],V_precond2[i+10000]);
		printf("at %d x is %d x_compare is  %f\n",i, x[i], x_compare[i]);
	}
	free(I);
	freee(J);
	free(V);
	free(I_precond);
	free(J_precond);
	free(V_precond);
	free(I_precondP);
	free(J_precondP);
	free(V_precondP);
	free(lowerJ);
	free(lowerI);
	free(lowerV);
	free(x);
	free(y);
	free(y_rodr);
	
	free(numInRow);
	free(numInRowL);
	free(numInRowLP);
	free(row_idx);
	free(row_idxL);
	free(row_idxLP);
	free(x_compare);
	free(diag);
	//interval2=(end_time2-start_time2)*1000/CLOCKS_PER_SEC;

	//printf("time consuming CPU is %f, time consuming GPU is %f, speedup is %f\n", interval1, interval2, interval1/interval2);
	//float Gflop=(totalNum*4+12*dimension)/interval1*1000*MAXIter;
	//printf("error is %f, total num is %d, time is %f ms, Gflops is %f, final error is %f\n",result_error/dimension, totalNum, interval1, Gflop, error_track[MAXIter-1]*1000);
	return 0;
}






