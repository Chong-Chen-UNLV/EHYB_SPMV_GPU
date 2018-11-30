#include "fspai.h"
#include "mmio.h"
#include "kernel.h"
#include "reordering.h"
#include "solver.h"
#include <unistd.h>


int main(int argc, char* argv[])
{
	unsigned int MAXthread = 2;	
	unsigned int totalNum=0;
	unsigned int lowerNum=0;
	unsigned int lowerNumPrecond=0;
	unsigned int totalNumPrecond=0;
	unsigned int totalNumPrecondP=0;
	int ret_code;
	MM_typecode matcode,matcode2;
	unsigned int MAXIter = 0;
	FILE *f,*f2;

	unsigned int dimension, N, bandwidth,bandwidthPrecond;   
	unsigned int maxRowNum, maxRowNumL, maxRowNumLP;
	unsigned int *I, *J;
	unsigned int *lowerI, *lowerJ;
	double *V;
	double *lowerV;
	double *x;
	double *y;
	double *x_compare;
	double *error_track;
	double result_error;
	int finish;
	char fileName[100];
	fileName[0] = '\0';
	char fileName2[100];
	unsigned int blocks;
	int oc;
	int arg_val;
	char* arg_str;
	cb_s cb;
    init_cb(&cb);
			
	while ((oc = getopt(argc, argv, "m:c:r:g:t:")) != -1) {
		switch (oc) {
			case 'm':
				/* input matrix */
				sprintf(fileName, "../read/%s.mtx", optarg);		
				printf("filename is %s\n", fileName);
				break;
			case 'c':
				/* the number of cycles */
				MAXIter = atoi(optarg);
				break;
			case 't':
				/* the number of cycles */
				MAXthread = atoi(optarg);
				break;
			case 'r':
				if(atoi(optarg) == 1)
					cb.RODR = true;
				break;
			case 'b':
				if(atoi(optarg) == 1)
					cb.BLOCK = true;
			case 'g':
				if(atoi(optarg) == 1)
					cb.GPU = true;
				break;
			case ':':
				       /* error handling, see text */
				printf("missing arguments\n");
				exit(0);
				break;
			case '?':
				printf("unrecongnized option\n");
				break;
			default:
				printf("option/arguments error!\n");       /* error handling, see text */
				exit(0);
		}
	}
	if (fileName[0] == '\0' || MAXIter == 0){
		printf("file name or max iteration number missing\n");
		exit(0);
	}

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
	int _dimension, _N, _lowerNum;
	if ((ret_code = mm_read_mtx_crd_size(f, &_dimension, &_N, &_lowerNum)) !=0)
		exit(1);	
	dimension = _dimension,
	N = _N;
	lowerNum = _lowerNum;

	/*The overall number of nozeros in this matrix*/
	totalNum=lowerNum*2-dimension;
	totalNumPrecond=lowerNum;
	totalNumPrecondP = lowerNum;

	size_t size0=lowerNum*sizeof(unsigned int);
	size_t size1=lowerNum*sizeof(double);
	size_t size2=totalNum*sizeof(unsigned int);
	size_t size3=totalNum*sizeof(double);
	size_t size4=dimension*sizeof(unsigned int);
	size_t size5=dimension*sizeof(double);
	size_t size6=lowerNum*sizeof(unsigned int);
	size_t size7=lowerNum*sizeof(double);


	lowerJ=(unsigned int *) malloc(size0);
	lowerI=(unsigned int *) malloc(size0);
	lowerV=(double *) malloc(size1);
	I=(unsigned int *) malloc(size2);
	J=(unsigned int *) malloc(size2);
	V=(double *) malloc(size3);
	x=(double *) malloc(size5);
	y=(double *) malloc(size5);
	double *diag=(double *) malloc(size5);
	x_compare=(double *) malloc(size5);
	error_track=(double *) malloc(MAXIter*sizeof(double));

	unsigned int *numInRowL;
	unsigned int *row_idxL;
	unsigned int *numInRowLP;
	unsigned int *row_idxLP;	
	unsigned int *numInRow;

	numInRowL=(unsigned int *) malloc(size4);
	numInRowLP=(unsigned int *) malloc(size4);
	row_idxL=(unsigned int *) malloc(size4 + sizeof(unsigned int));
	row_idxLP=(unsigned int *) malloc(size4 + sizeof(unsigned int));
	numInRow=(unsigned int *) malloc(size4);
	unsigned int tempI, tempJ;
	double tempV;
	for (int i=0; i<lowerNum; i++)
	{
		fscanf(f, "%d %d %lg\n", &tempI, &tempJ, &tempV);
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


	unsigned int *row_idx=(unsigned int *)malloc((dimension+1)*sizeof(unsigned int));
	maxRowNum=0;
	maxRowNumL=0;
	maxRowNumLP=0;
	row_idx[0] = 0;
	row_idxL[0] = 0;
	row_idxLP[0] = 0;
	for (int i=1;i<= dimension;i++)
	{

		if (numInRow[i-1]>maxRowNum)
			maxRowNum=numInRow[i-1];
		if (numInRowL[i-1]>maxRowNumL)
			maxRowNumL=numInRowL[i-1];
		if (numInRowLP[i-1]>maxRowNumLP)
			maxRowNumLP=numInRowLP[i-1];			

		row_idx[i]=row_idx[i-1]+numInRow[i-1];
		row_idxL[i]=row_idxL[i-1]+numInRowL[i-1];
		row_idxLP[i]=row_idxLP[i-1]+numInRowLP[i-1];
		numInRow[i-1]=0;
		numInRowLP[i-1]=0;
		//determine y
	}	
	if (numInRow[dimension-1]>maxRowNum) maxRowNum=numInRow[dimension-1];
	if (numInRowL[dimension-1]>maxRowNumL) maxRowNumL=numInRowL[dimension-1];	
	if (numInRowLP[dimension-1]>maxRowNumLP) maxRowNumLP=numInRowLP[dimension-1];
	printf("maxRowNum is %d, maxRowNumPr is %d, maxRowNumLP is %d\n", maxRowNum, maxRowNumL, maxRowNumLP);
	numInRow[dimension-1]=0;
	numInRowLP[dimension-1]=0;
	for (int i=0;i<dimension;i++)
	{		
		srand(i);
		//x_compare[i]=(double) (rand()%200-100)/100;
		x_compare[i]=1;
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
	unsigned int *I_rodr, *J_rodr, *part_boundary, *rodr_list;
	double *V_rodr, *x_rodr, *y_rodr;
	
	if(cb.RODR){
		blocks = ceil(dimension/(shared_per_block/element_size));
		printf("blocks is %d\n", blocks);
		rodr_list = (unsigned int* )calloc(dimension, sizeof(unsigned int)); 
		part_boundary = (unsigned int* )calloc((blocks + 1), sizeof(unsigned int)); 	
		I_rodr = (unsigned int *) malloc(size2);
		J_rodr = (unsigned int *) malloc(size2);
		V_rodr = (double *) malloc(size3);
		x_rodr = (double* )calloc(dimension, sizeof(double)); 
		y_rodr = (double* )calloc(dimension, sizeof(double)); 
		/*NOTICE: maxRowNum and maxRowNum for preconditioners do not needed 
		to be updated. When we do reordering, all elements of certain row
		will be assigned to SAME row according to rodr_list
		so the maxRowNum should be same, only occured in different row number
		*/
		matrix_reorder(&dimension, totalNum, I, J, V, numInRow, row_idx,
				 I_rodr, J_rodr, V_rodr, rodr_list, part_boundary, blocks);
		vector_reorder(dimension, y, y_rodr, rodr_list);
		update_numInRowL(totalNum, 
			dimension, 
			I_rodr, 
			J_rodr, 
			V_rodr, 
			numInRowL,
			row_idxL, 
			row_idxLP,
			diag);
	}
	/*---------------------read the preconditioner ------------------------------*/

	unsigned int *I_precond=(unsigned int *) malloc(size6);
	unsigned int *J_precond=(unsigned int *) malloc(size6);
	double* V_precond=(double*) malloc(size7);	

	/*int rt=pthread_barrier_init(&barr, NULL, MAXthread);
	  rt=pthread_barrier_init(&barr1, NULL, MAXthread);*/

	S Sthread[MAXthread];

	for (int t=0;t<MAXthread;t++)
	{
		if(cb.RODR){
			Sthread[t].I=I_rodr;
			Sthread[t].J=J_rodr;
			Sthread[t].V= V_rodr;
		}
		else{
			Sthread[t].I=I;
			Sthread[t].J=J;
			Sthread[t].V= V;
		}
		Sthread[t].I_precond=I_precond;
		Sthread[t].J_precond=J_precond;
		Sthread[t].V_precond=V_precond;
		Sthread[t].maxRowNum=maxRowNum;
		Sthread[t].numInRow=numInRow;
		Sthread[t].row_idx=row_idx;
		Sthread[t].numInRowPrecond=numInRowL;
		Sthread[t].row_idxPrecond=row_idxL;
		Sthread[t].diag = diag;
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
	unsigned int *I_precondP=(unsigned int *) malloc(size6);
	unsigned int *J_precondP=(unsigned int *) malloc(size6);
	double *V_precondP=(double *) malloc(size7);
	
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
	
	unsigned int realIter, pp, procNum;
	struct timeval start, end;

	#pragma omp parallel
	{
		#pragma omp master
		{
			procNum=omp_get_num_threads();
		}
	}
	printf("thread num is %d\n",procNum);
	if(cb.GPU){
        /*please notice that we need transfer the format of matrix to HYB, so we need I, J, V completely
	   	for GPU solver, for CPU solver, no format change is applied, so we only need row_idx**/
        matrixCOO_S matrix, matrix_precond, matrix_precondP;
        if(cb.RODR){
            init_matrixCOO_S(&matrix, dimension, totalNum, maxRowNum, row_idx, numInRow, I_rodr, J_rodr, V_rodr);
        }
        else{
            int* part_boundary = NULL;
            init_matrixCOO_S(&matrix, dimension, totalNum, maxRowNum, row_idx, numInRow, I, J, V);
        }
        init_matrixCOO_S(&matrix_precond, dimension, totalNumPrecond, maxRowNumL, row_idxL, numInRowL, I_precond, J_precond, V_precond);
        init_matrixCOO_S(&matrix_precondP, dimension, totalNumPrecondP, maxRowNumLP, row_idxLP, numInRowLP, I_precondP, J_precondP, V_precondP);
        if(cb.RODR){
            solverGPU_HYB(&matrix, &matrix_precond, &matrix_precondP,
                    y_rodr, x_rodr, MAXIter, &realIter, cb, blocks, part_boundary);
		}
		else{
            solverGPU_HYB(&matrix, &matrix_precond, &matrix_precondP,
                    y_rodr, x_rodr, MAXIter, &realIter, cb, blocks, part_boundary);
		}
	}
	else if(!(cb.GPU)){
		if(cb.RODR){
			solverPrecondCPU(procNum, dimension, totalNum, row_idx, J_rodr, V_rodr, totalNumPrecond, 
					row_idxL, J_precond, V_precond,totalNumPrecondP, 
					row_idxLP, J_precondP, V_precondP, 
					y_rodr, x_rodr, MAXIter, &realIter);
		}
		else if(!cb.RODR){
			solverPrecondCPU(procNum, dimension, totalNum, row_idx, J, V, totalNumPrecond, 
					row_idxL, J_precond, V_precond,totalNumPrecondP, 
					row_idxLP, J_precondP, V_precondP, 
					y, x, MAXIter, &realIter);
		}
		else;
	}
	else;
	//solverPrecondGPU_CUSPARSE(dimension, totalNum, *row_idx, *J, *V, totalNumPrecond, *row_idxL, *J_precond, *V_precond, 
	//	totalNumPrecondP, *row_idxLP, *J_precondP, *V_precondP, *y, *x, MAXIter, *realIter);
		
	if(cb.RODR){
		vector_recover(dimension,  x_rodr, x, rodr_list);
	}

	for (int i=0;i<10;i++)
	{
		//printf("Xeon_phi I is %d J %d is V is %f\n",I_precond[i+10000], J_precond[i+10000], V_precond[i+10000]);
		//printf("CPU I is %d, J is %d, V is %f\n",I_precond2[i+10000],J_precond2[i+10000],V_precond2[i+10000]);
		printf("at %d x is %f x_compare is  %f\n",i, x[i], x_compare[i]);
	}
	free(I);
	free(J);
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
	if(cb.RODR){
		free(y_rodr);
		free(I_rodr);
		free(J_rodr);
		free(V_rodr);
	}
	
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
	//double Gflop=(totalNum*4+12*dimension)/interval1*1000*MAXIter;
	//printf("error is %f, total num is %d, time is %f ms, Gflops is %f, final error is %f\n",result_error/dimension, totalNum, interval1, Gflop, error_track[MAXIter-1]*1000);
	return 0;
}






