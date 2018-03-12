#include "fspai.h"
#include "mmio.h"
#include <omp.h>

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
	int *rowNumAccumL;
	int *numInRowLP;
	int *rowNumAccumLP;	
	int *numInRow;

	numInRowL=(int *) malloc(size4);
	numInRowLP=(int *) malloc(size4);
	rowNumAccumL=(int *) malloc(size4 + sizeof(int));
	rowNumAccumLP=(int *) malloc(size4 + sizeof(int));
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


	int *rowNumAccum=(int *)malloc((dimension+1)*sizeof(int));
	maxRowNum=0;
	maxRowNumPrecond=0;
	maxRowNumPrecondP=0;
	rowNumAccum[0] = 0;
	rowNumAccumL[0] = 0;
	rowNumAccumLP[0] = 0;
	for (int i=1;i<= dimension;i++)
	{

		if (numInRow[i-1]>maxRowNum)
			maxRowNum=numInRow[i-1];
		if (numInRowL[i-1]>maxRowNumPrecond)
			maxRowNumPrecond=numInRowL[i-1];
		if (numInRowLP[i-1]>maxRowNumPrecondP)
			maxRowNumPrecondP=numInRowLP[i-1];			

		rowNumAccum[i]=rowNumAccum[i-1]+numInRow[i-1];
		rowNumAccumL[i]=rowNumAccumL[i-1]+numInRowL[i-1];
		rowNumAccumLP[i]=rowNumAccumLP[i-1]+numInRowLP[i-1];
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
		//x_compare[i]=(float) (rand()%200-100)/100;
		x_compare[i]=1;
	}
	int index1, index2;

	for (int i=0;i<lowerNum;i++)
	{
		tempI=lowerI[i];
		tempJ=lowerJ[i];
		tempV=lowerV[i];
		index1=rowNumAccum[tempI]+numInRow[tempI];
		index2=rowNumAccum[tempJ]+numInRow[tempJ];
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
		



	/*---------------------read the preconditioner ------------------------------*/


	int *I_precond=(int *) malloc(size6);
	int *J_precond=(int *) malloc(size6);
	float *V_precond=(float *) malloc(size7);	

	/*int rt=pthread_barrier_init(&barr, NULL, MAXthread);
	  rt=pthread_barrier_init(&barr1, NULL, MAXthread);*/

	S Sthread[MAXthread];

	for (int t=0;t<MAXthread;t++)
	{
		Sthread[t].I=I;
		Sthread[t].J=J;
		Sthread[t].V=V;
		Sthread[t].I_precond=I_precond;
		Sthread[t].J_precond=J_precond;
		Sthread[t].V_precond=V_precond;
		Sthread[t].maxRowNum=maxRowNum;
		Sthread[t].numInRow=numInRow;
		Sthread[t].rowNumAccum=rowNumAccum;
		Sthread[t].numInRowPrecond=numInRowL;
		Sthread[t].rowNumAccumPrecond=rowNumAccumL;
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
		if (tempI<0||tempI>dimension-1||tempJ<0||tempJ>dimension-1) printf("error happend at %d with tempI %d and tempJ %d\n", i ,tempI, tempJ); 
		index1=rowNumAccumLP[tempI]+numInRowLP[tempI];
		I_precondP[index1]=tempI;
		J_precondP[index1]=tempJ;
		V_precondP[index1]=tempV;
		numInRowLP[tempI]+=1;
	}


	int *I_accum, *I_precond_accum, *I_precondP_accum;
	I_accum = (int *)malloc((dimension + 1)*sizeof(int));
	I_precond_accum = (int *)malloc((dimension + 1)*sizeof(int));
	I_precondP_accum = (int *)malloc((dimension + 1)*sizeof(int));
	for (int i = 0; i <= dimension; ++i){
		I_accum[i] = rowNumAccum[i];
		I_precond_accum[i] = rowNumAccumL[i];
		I_precondP_accum[i] = rowNumAccumLP[i];
	}

	int realIter;



	struct timeval start, end;
	int pp;


	

	int procNum;

#pragma omp parallel
	{
#pragma omp master
		{
			procNum=omp_get_num_threads();
		}
	}
	printf("thread num is %d\n",procNum);
	//solverPrecondCPU(dimension, totalNum, *I_accum, *J, *V, totalNumPrecond, *I_precond_accum, *J_precond, *V_precond, 
	//	totalNumPrecondP, *I_precondP_accum, *J_precondP, *V_precondP, *y, *x, MAXIter, *realIter);
	solverPrecondGPU_Chong(dimension, totalNum, *I_accum, *J, *V, totalNumPrecond, *I_precond_accum, *J_precond, *V_precond, 
		totalNumPrecondP, *I_precondP_accum, *J_precondP, *V_precondP, *y, *x, MAXIter, *realIter);
	//solverPrecondGPU_CUSPARSE(dimension, totalNum, *I_accum, *J, *V, totalNumPrecond, *I_precond_accum, *J_precond, *V_precond, 
	//	totalNumPrecondP, *I_precondP_accum, *J_precondP, *V_precondP, *y, *x, MAXIter, *realIter);
	

	for (int i=0;i<10;i++)
	{
		//int base=rowNumAccumL[I[a]];
		//printf("Xeon_phi I is %d J %d is V is %f\n",I_precond[i+10000], J_precond[i+10000], V_precond[i+10000]);
		//printf("CPU I is %d, J is %d, V is %f\n",I_precond2[i+10000],J_precond2[i+10000],V_precond2[i+10000]);
		//printf("Y %d is %f\n",i,y[i]);
	}


	//interval2=(end_time2-start_time2)*1000/CLOCKS_PER_SEC;

	//printf("time consuming CPU is %f, time consuming GPU is %f, speedup is %f\n", interval1, interval2, interval1/interval2);
	//float Gflop=(totalNum*4+12*dimension)/interval1*1000*MAXIter;
	//printf("error is %f, total num is %d, time is %f ms, Gflops is %f, final error is %f\n",result_error/dimension, totalNum, interval1, Gflop, error_track[MAXIter-1]*1000);
	return 0;
}






