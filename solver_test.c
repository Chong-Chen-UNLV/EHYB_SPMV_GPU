#include "kernel.h"
#include "reordering.h"
#include "spmv.h"
#include <unistd.h>
#include "mmio.h"


static void compare(float* yResult, float* y, const float threshold, const int dimension)
{
	float avgdiff = 0; 
	float avgampldiff = 0;
	int k = 0;
	for (int i = 0; i < dimension; ++i)
	{
		float d = fabs(y[i] - yResult[i]);
		float ampl = fmin(fabs(y[i]), fabs(yResult[i]));
		if (d > ampl*threshold)
		{
			if(k < 100){
				printf("large difference at %d  : realy %f vs yResult %f\n", i, y[i] , yResult[i]);
				k++;
			}
		}
		avgdiff += d;
		if(ampl > 0) avgampldiff += d/ampl;
	}
	avgdiff ;
	avgampldiff ;
	printf("diff is %e, ampldiff is %e\n", avgdiff, avgampldiff);
}

static int matrixRead_unsym(matrixCOO* localMatrixCOO, float** xCompare_in, float** y_in, FILE *f)
{ 	int _dimension, _N;
	int ret_code, totalNum;
	if ((ret_code = mm_read_mtx_crd_size(f, &_dimension, &_N, &totalNum)) !=0)
		exit(1);	
	localMatrixCOO->dimension = _dimension;
	*xCompare_in = (float*)malloc(_dimension*sizeof(float));
	*y_in = (float*)malloc(_dimension*sizeof(float));
	float*	xCompare = *xCompare_in;
	float*	y = *y_in;
	localMatrixCOO->totalNum = totalNum;
	localMatrixCOO->nParts = ceil(((float) _dimension)/vectorCacheSize);
	printf("parts is %d with cachSize %d\n", localMatrixCOO->nParts, vectorCacheSize);
	if(localMatrixCOO->nParts <= 40){
		localMatrixCOO->kernelPerPart = smSize/(localMatrixCOO->nParts);
		printf("kernel per part is %d\n", localMatrixCOO->kernelPerPart);
	}
	localMatrixCOO->partBoundary = (int* )calloc(_dimension, sizeof(int));
	localMatrixCOO->reorderList = (int* )calloc(_dimension, sizeof(int));
	localMatrixCOO->numInRow = (int* )calloc(_dimension, sizeof(int));
	localMatrixCOO->numInRow2 = (int* )calloc(_dimension, sizeof(int));

	localMatrixCOO->I=(int *) malloc(totalNum*sizeof(int));
	localMatrixCOO->J=(int *) malloc(totalNum*sizeof(int));
	localMatrixCOO->V=(float *) malloc(totalNum*sizeof(float));
	localMatrixCOO->diag = (float *) malloc(_dimension*sizeof(float));

	int* numInRow = localMatrixCOO->numInRow;
	int* I = localMatrixCOO->I;
	int* J = localMatrixCOO->J;
	float* V = localMatrixCOO->V;

	for (int i=0;i < _dimension;i++){
		srand(i);
		xCompare[i]=(float) (rand()%200-100)/1000;
		//x_compare[i]=1;
	}
	int tempI, tempJ;
	float tempV;
	for (int i=0; i<totalNum; i++){
		fscanf(f, "%d %d %g\n", &tempI, &tempJ, &tempV);
		J[i]=tempJ-1;  /* adjust from 1-based to 0-based */
		I[i]=tempI-1;
		V[i]=tempV;
		numInRow[tempI-1]+=1;
		y[tempI-1]+=tempV*xCompare[tempJ-1];
	}

	localMatrixCOO->rowIdx=(int *)malloc((_dimension+1)*sizeof(int));
	int* rowIdx = localMatrixCOO->rowIdx;
	int maxCol = 0;
	int maxColIdx = 0;
	rowIdx[0] = 0;
	
	for (int i=1;i<= _dimension;i++){
		if (numInRow[i-1]>maxCol){
			maxCol = numInRow[i-1];
			maxColIdx = i-1;
		}

		rowIdx[i]=rowIdx[i-1]+numInRow[i-1];
	}	
	if (numInRow[_dimension-1] > maxCol) maxCol=numInRow[_dimension-1];
	printf("maxCol is %d at idx %d\n", maxCol, maxColIdx);
	localMatrixCOO->maxCol = maxCol;
	return 1;
}
static int matrixRead_sym(matrixCOO* localMatrixCOO, float** xCompare_in, float** y_in, FILE *f)
{
	int _dimension, _N, _lowerNum;
	int ret_code, totalNum, lowerNum;
	if ((ret_code = mm_read_mtx_crd_size(f, &_dimension, &_N, &_lowerNum)) !=0)
		exit(1);	
	localMatrixCOO->dimension = _dimension;
	lowerNum = _lowerNum;
	totalNum = lowerNum*2-_dimension;
	
	*xCompare_in = (float*)malloc(_dimension*sizeof(float));
	*y_in = (float*)malloc(_dimension*sizeof(float));
	float*	xCompare = *xCompare_in;
	float*	y = *y_in;
	/*The overall number of nozeros in this matrix*/
	localMatrixCOO->totalNum = totalNum;
	localMatrixCOO->nParts = ceil(((float) _dimension)/vectorCacheSize);
	printf("parts is %d with cahe size %d\n", localMatrixCOO->nParts, vectorCacheSize);
	if(localMatrixCOO->nParts <= 40){
		localMatrixCOO->kernelPerPart = smSize/(localMatrixCOO->nParts);
		printf("kernel per part is %d\n", localMatrixCOO->kernelPerPart);
	}
	localMatrixCOO->partBoundary = (int* )calloc(_dimension, sizeof(int));
	localMatrixCOO->reorderList = (int* )calloc(_dimension, sizeof(int));
	localMatrixCOO->numInRow = (int* )calloc(_dimension, sizeof(int));
	localMatrixCOO->numInRow2 = (int* )calloc(_dimension, sizeof(int));

	int* lowerI=(int *) malloc(lowerNum*sizeof(int));
	int* lowerJ=(int *) malloc(lowerNum*sizeof(int));
	float* lowerV=(float *) malloc(lowerNum*sizeof(float));

	localMatrixCOO->I=(int *) malloc(totalNum*sizeof(int));
	localMatrixCOO->J=(int *) malloc(totalNum*sizeof(int));
	localMatrixCOO->V=(float *) malloc(totalNum*sizeof(float));
	localMatrixCOO->diag = (float *) malloc(_dimension*sizeof(float));
	int* numInRow = localMatrixCOO->numInRow;
	int* I = localMatrixCOO->I;
	int* J = localMatrixCOO->J;
	float* V = localMatrixCOO->V;
	float* diag = localMatrixCOO->diag;
	
	int tempI, tempJ;
	float tempV;
	for (int i=0; i<lowerNum; i++){
		fscanf(f, "%d %d %f\n", &tempI, &tempJ, &tempV);
		lowerJ[i]=tempJ-1;  /* adjust from 1-based to 0-based */
		lowerI[i]=tempI-1;
		lowerV[i]=tempV;
		numInRow[tempI-1]+=1;
		if (tempI!=tempJ)
		{
			numInRow[tempJ-1]+=1;
		}		
	}

	localMatrixCOO->rowIdx=(int *)malloc((_dimension+1)*sizeof(int));
	int* rowIdx = localMatrixCOO->rowIdx;
	int maxCol = 0;
	int maxColIdx = 0;
	rowIdx[0] = 0;
	
	for (int i=1;i<= _dimension;i++){
		if (numInRow[i-1]>maxCol){
			maxCol = numInRow[i-1];
			maxColIdx = i-1;
		}

		rowIdx[i]=rowIdx[i-1]+numInRow[i-1];
		numInRow[i-1]=0;
	}	
	if (numInRow[_dimension-1] > maxCol) maxCol=numInRow[_dimension-1];
	
	printf("maxCol is %d with idx %d\n", maxCol, maxColIdx);
	localMatrixCOO->maxCol = maxCol;
	numInRow[_dimension-1]=0;
	for (int i=0;i < _dimension;i++){
		srand(i);
		xCompare[i]=(float) (rand()%200-100)/1000;
		//x_compare[i]=1;
	}
	int index1, index2;

	for (int i=0;i<lowerNum;i++){
		tempI=lowerI[i];
		tempJ=lowerJ[i];
		tempV=lowerV[i];
		if(tempJ >= _dimension || tempI >= _dimension)
			exit(0);
		index1=rowIdx[tempI]+numInRow[tempI];
		index2=rowIdx[tempJ]+numInRow[tempJ];
		numInRow[tempI]+=1;
		I[index1]=tempI;
		J[index1]=tempJ;
		V[index1]=tempV;
		y[tempI]+=tempV*xCompare[tempJ];
		if (tempI != tempJ)
		{
			numInRow[tempJ]+=1;
			I[index2]=tempJ;
			J[index2]=tempI;
			V[index2]=tempV;
			y[tempJ]+=tempV*xCompare[tempI];
		}
		else
		{
			diag[tempI]=tempV;
		}
	}
	free(lowerV);
	free(lowerI);
	free(lowerJ);
	return 1;
}

int main(int argc, char* argv[])
{
	static MM_typecode matcode;
	int MAXIter = 0;
	FILE *f;
	float *y;
	float *xCompare;
	char fileName[100];
	fileName[0] = '\0';
	int oc;
	cb_s cb;
    init_cb(&cb);
			
	while ((oc = getopt(argc, argv, "m:i:r:t:f:p:")) != -1) {
		switch (oc) {
			case 'm':
				/* input matrix */
				sprintf(fileName, "../read/%s.mtx", optarg);		
				printf("filename is %s\n", fileName);
				break;
			case 'i':
				/* the number of cycles */
				MAXIter = atoi(optarg);
				break;
			case 't':
				/* the number of threads*/
				//MAXthread = atoi(optarg);
				break;
			case 'p':
				if(atoi(optarg) == 1)
					cb.PRECOND = true;
				break;
		
			case 'f':
				/*using SPAI instaed of factorized SPAI*/
				if(atoi(optarg) == 1)
					cb.FACT = false;
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
	if (!cb.RODR || !cb.CACHE || !cb.BLOCK){
		printf("this program only test RODR, BLOCK, and CACHE enabled case\n");
		exit(0);
	}

	//---------------------------------read the matrix---------------------------
	if ((f = fopen(fileName, "r")) == NULL){ 
		printf("file read error\n");
		exit(1);
	}

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

	matrixCOO localMatrixCOO;
	if (mm_is_symmetric(matcode)){
		printf("read symmetric matrix\n");
		matrixRead_sym(&localMatrixCOO, &xCompare, &y, f);
	} else {
		printf("read unsymmetric matrix\n");
		matrixRead_unsym(&localMatrixCOO, &xCompare, &y, f);
	}
	fclose(f);
	float *yResult = (float *) calloc(localMatrixCOO.dimension, sizeof(float));

	//spmvHYB(&localMatrixCOO, xCompare, yResult, MAXIter);
	//spmvGeneric(&localMatrixCOO, xCompare, yResult, MAXIter);
	//solverGPuUnprecondCUSPARSE(&localMatrixCOO, xCompare, yResult, MAXIter);	
	//for (int i=0;i<10;i++)
	//{
	//	printf("at %d yResult is %f y is  %f\n",i + 30000, yResult[i + 30000], y[i + 30000]);
	//}
	//return 0;

	float *xReorder = (float* )calloc(localMatrixCOO.dimension, sizeof(float)); 
	float *yReorder = (float* )calloc(localMatrixCOO.dimension, sizeof(float)); 
	if (mm_is_symmetric(matcode)){
		matrixReorder(&localMatrixCOO);
	} else {
		printf("unsymmetric reordering\n");
		matrixReorder_unsym(&localMatrixCOO);
	}

	vectorReorder(localMatrixCOO.dimension, xCompare, xReorder, localMatrixCOO.reorderList);

	int realIter; 

	//format change is completed in the solver function

	spmvGPuEHYB(&localMatrixCOO, xReorder, yReorder, MAXIter, &realIter);
	vectorRecover(localMatrixCOO.dimension, yReorder, yResult, localMatrixCOO.reorderList);

	for (int i=0;i<10;i++)
	{
		printf("at %d yResult is %f y is  %f\n",i + 30000, yResult[i + 30000], y[i + 30000]);
	}
	compare(yResult, y, 0.01, localMatrixCOO.dimension);
	free(localMatrixCOO.I);
	free(localMatrixCOO.J);
	free(localMatrixCOO.V);
	free(yResult);
	free(xReorder);
	free(y);
	
	free(localMatrixCOO.numInRow);
	free(localMatrixCOO.numInRow2);
	free(localMatrixCOO.rowIdx);
	free(xCompare);
	free(localMatrixCOO.diag);
	//interval2=(end_time2-start_time2)*1000/CLOCKS_PER_SEC;

	//printf("time consuming CPU is %f, time consuming GPU is %f, speedup is %f\n", interval1, interval2, interval1/interval2);
	//float Gflop=(totalNum*4+12*dimension)/interval1*1000*MAXIter;
	//printf("error is %f, total num is %d, time is %f ms, Gflops is %f, final error is %f\n",result_error/dimension, totalNum, interval1, Gflop, error_track[MAXIter-1]*1000);
	return 0;
}
