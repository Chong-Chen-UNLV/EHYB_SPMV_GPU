#include "fspai.h"
#include "mmio.h"
#include <omp.h>

#define MAXthread 16 
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
	rowNumAccumL=(int *) malloc(size4);
	rowNumAccumLP=(int *) malloc(size4);
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


	int *rowNumAccum=(int *)malloc(dimension*sizeof(int));
	maxRowNum=0;
	maxRowNumPrecond=0;
	maxRowNumPrecondP=0;
	for (int i=1;i<dimension;i++)
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


	int *I_accum, *J_1;
	float *V_1;
	int *I_precond_accum, *J_precond_1;
	float *V_precond_1;
	int *I_precondP_accum, *J_precondP_1;
	float *V_precondP_1;

	int totalNum_1,totalNumPrecond_1,totalNumPrecondP_1;

	int realIter;

	formatChange(dimension, numInRow, &totalNum_1, I, J, V, &I_accum, &J_1, &V_1);
	formatChange(dimension, numInRowL, &totalNumPrecond_1, I_precond, J_precond, V_precond, &I_precond_accum, &J_precond_1, &V_precond_1);
	formatChange(dimension, numInRowLP, &totalNumPrecondP_1, I_precondP, J_precondP, V_precondP, &I_precondP_accum, &J_precondP_1, &V_precondP_1);


	struct timeval start, end;
	int pp;


	float *rk=(float *)malloc(size5);
	float *zk=(float *)malloc(size5);
	float *zk1=(float *)malloc(size5);
	float *pk=(float *)malloc(size5);
	float *bp=(float *)malloc(size5);

	float alphak = 0,betak = 0,dotp0 = 0, dotrz0 = 0,dotrz1 = 0,doth = 0,alphak_1 = 0,dotr0 = 0;

	//float dotz0,dotz0_compare;

	float error=10000;
	float threshold;

	for (int i=0;i<dimension;i++){
		zk[i]=0;
		zk1[i]=0;
		rk[i]=y[i];
	}

	int procNum;

#pragma omp parallel
	{
#pragma omp master
		{
			procNum=omp_get_num_threads();
		}
	}
	printf("thread num is %d\n",procNum);
	int *boundary=(int *)malloc((procNum+1)*sizeof(int));
	int *precond_boundary=(int *)malloc((procNum+1)*sizeof(int));
	int *precondP_boundary=(int *)malloc((procNum+1)*sizeof(int));
	boundary[0]=0;
	precond_boundary[0]=0;
	precondP_boundary[0]=0;
	boundary[procNum]=dimension;
	precond_boundary[procNum]=dimension;
	precondP_boundary[procNum]=dimension;

	int stride, stridePrecond,stridePrecondP;

	stride=ceil((float)totalNum_1/procNum);
	stridePrecond=ceil((float)totalNumPrecond_1/procNum);
	stridePrecondP=ceil((float)totalNumPrecondP_1/procNum);
	int bias, biasPrecond,biasPrecondP;

	bias=1;
	biasPrecond=1;
	biasPrecondP=1;

	int a,b;

	for (int i=0;i<dimension;i++){

		if (I_accum[i]>bias*stride){
			boundary[bias]=i;
			bias++;
		}

		if (I_precond_accum[i]>biasPrecond*stridePrecond){
			precond_boundary[biasPrecond]=i;
			biasPrecond++;
		}		

		if (I_precondP_accum[i]>biasPrecondP*stridePrecondP){
			precondP_boundary[biasPrecondP]=i;
			biasPrecondP++;
		}

	}

	__assume_aligned(V_1,64);
	__assume_aligned(J_1,64);;

	__assume_aligned(V_precond_1,64);
	__assume_aligned(J_precond_1,64);

	__assume_aligned(V_precondP_1,64);
	__assume_aligned(J_precondP_1,64);

	/*#pragma omp parallel for
	  for (int i=0;i<totalNumPrecond_1;i++)
	  zk1[I_precond_1[i]]+=V_precond_1[i]*rk[J_precond_1[i]];*/
#pragma omp parallel private(rank,tempV)
	{
		rank=omp_get_thread_num();
		for (int i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
			tempV=0;
			a=I_precond_accum[i];
			b=I_precond_accum[i+1];
#pragma simd
			for (int j=a;j<b;j++){
				tempV+=V_precond_1[j]*rk[J_precond_1[j]];
			}
			zk1[i]+=tempV;
			//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
		}
	}		

#pragma omp parallel private(rank,tempV)
	{
		rank=omp_get_thread_num();
		for (int i=precondP_boundary[rank];i<precondP_boundary[rank+1];i++){
			tempV=0;
			a=I_precondP_accum[i];
			b=I_precondP_accum[i+1];
#pragma simd
			for (int j=a;j<b;j++)
				tempV+=V_precondP_1[j]*zk1[J_precondP_1[j]];
			zk[i]+=tempV;
		}

	}				

	//initialize_all(dimension,pk_d,bp_d,x_d,zk_d,vector_in_d);
#pragma omp parallel for
	for (int i=0;i<dimension;i++){
		pk[i]=zk[i];
		bp[i]=0;
		x[i]=0;
	}


	//cublasSdot(handle, dimension,vector_in_d,1,vector_in_d,1,&doth);
	//#pragma omp parallel for reduction(+:doth)
	for (int i=0;i<dimension;i++){
		doth+=y[i]*y[i];
	}
	float errorNorm=sqrt(doth);
	//printf("errorNorm is %f\n",errorNorm);
	int iter=0;


	gettimeofday(&start1, NULL);
	/*-----------Start the iteration by a while loop-----------*/
	while (iter<MAXIter&&error>1E-7){
		//matrix-vector multiplication, accelerated by our own functions
		/*#pragma omp parallel for
		  for (int i=0;i<totalNum_1;i++)
		  bp[I_1[i]]+=V_1[i]*pk[J_1[i]];*/
#pragma omp parallel private(rank,tempV,start2, end2)
		{
			rank=omp_get_thread_num();
			//gettimeofday(&start2, NULL);
			for (int i=boundary[rank];i<boundary[rank+1];i++){
				tempV=0;
				a=I_accum[i];
				b=I_accum[i+1];
#pragma simd
				for (int j=a;j<b;j++){
					tempV+=V_1[j]*pk[J_1[j]];
				}
				bp[i]+=tempV;
			}
			//gettimeofday(&end2, NULL);
			//printf("at iter %d, rank is %d, boundary span is %d, time is %d us\n", iter, rank, boundary[rank+1]-boundary[rank], ((end2.tv_sec * 1000000 + end2.tv_usec)-(start2.tv_sec * 1000000 + start2.tv_usec)));
		}

		dotp0=0;
#pragma omp parallel for reduction(+:dotp0)
		for (int i=0;i<dimension;i++)
			dotp0+=bp[i]*pk[i];

		//r_k*z_k
		dotrz0=0;
#pragma omp parallel for reduction(+:dotrz0)
		for (int i=0;i<dimension;i++)
			dotrz0+=rk[i]*zk[i];

		alphak=dotrz0/dotp0;
		alphak_1=-alphak;

		//update x_k to x_k+1, x=x+alphak*p;
		//update r_k to r_k+1, r=r-alphak*bp;
#pragma omp parallel for
		for (int i=0;i<dimension;i++){
			x[i]+=alphak*pk[i];
			rk[i]+=alphak_1*bp[i];
			zk1[i]=0;
			zk[i]=0;
		}

		//SpMV: zk=inv(m)*rk
#pragma omp parallel private(rank,tempV)
		{
			rank=omp_get_thread_num();
			for (int i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
				tempV=0;
				a=I_precond_accum[i];
				b=I_precond_accum[i+1];
#pragma simd
				for (int j=a;j<b;j++){
					tempV+=V_precond_1[j]*rk[J_precond_1[j]];
				}
				zk1[i]+=tempV;
				//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
			}
		}		

#pragma omp parallel private(rank,tempV)
		{
			rank=omp_get_thread_num();
			for (int i=precondP_boundary[rank];i<precondP_boundary[rank+1];i++){
				tempV=0;
				a=I_precondP_accum[i];
				b=I_precondP_accum[i+1];
#pragma simd
				for (int j=a;j<b;j++)
					tempV+=V_precondP_1[j]*zk1[J_precondP_1[j]];
				zk[i]+=tempV;
			}

		}			

		//r_k+1 * z_k+1
		dotrz1=0;
#pragma omp parallel for reduction (+:dotrz1)
		for (int i=0;i<dimension;i++)
			dotrz1+=rk[i]*zk[i];
		betak=dotrz1/dotrz0;
		//printf("dotp0 is %f, dotrz0 is %f dotrz1 is %f, betak is %f and alphak is %f at iter %d\n", dotp0, dotrz0,dotrz1,betak, alphak, iter);
		//p=r+gamak*p;
#pragma omp parallel for
		for (int i=0;i<dimension;i++){
			pk[i]=zk[i]+betak*pk[i];
			bp[i]=0;
		}
		dotr0=0;
#pragma omp parallel for reduction (+:dotr0)
		for (int i=0;i<dimension;i++)
			dotr0+=rk[i]*rk[i];
		//printf("dotr0 is %f\n",dotrz1);		
		error=sqrt(dotr0)/errorNorm;
		/*#pragma omp master
		  printf("at iter %d, alpha is %f, beta is %f, dotrz0 is %f, dotrz1 is %f, dotp0 is %f\n", iter, alphak, betak, dotrz0, dotrz1, dotp0);*/
		iter++;
	}
	/*--------------finish the iterations---------------------*/
	gettimeofday(&end1, NULL);	
	printf("max iteration is %d with error %e\n",iter, error);

	/*------------------------------------------------------------Finish the iterative solver-------------------------------------------------*/



	gettimeofday(&end1, NULL);	
	printf("Solver Xeon_phi time is %ld us\n",(end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec));
	float timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, Xeon_Phi Gflops is %f\n ",iter, (1e-9*(totalNum*4+14*dimension)*1000*iter)/timeByMs);



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

void insertSort(int *J, float *V, int num, int *outJ, float *outV){
	//input is the unsorted J point, output is the sorted J and V

	int tempJ,tempMap;
	int i,j,k;
	float tempV;
	for (i=0;i<num;i++){
		outJ[i]=0;
		outV[i]=0;
	}
	//initial output[0]
	outJ[0]=*J;
	outV[0]=*V;
	for (i=1;i<num;i++){
		tempJ=J[i];
		tempV=V[i];
		//updata output [i]
		for (j=i;j>0;j--){
			if (tempJ>=outJ[j-1]){
				outV[j]=tempV;
				outJ[j]=tempJ;
				break;
			}
			else
			{
				outJ[j]=outJ[j-1];
				outV[j]=outV[j-1];
				//output[j-1]=*(input+i);
				//outputDegree[j-1]=tempDegree;
			}
		}
		if (j==0)
		{		
			outJ[0]=tempJ;
			outV[0]=tempV;
		}
	}
}

void solver(const int dimension, const int totalNum, const int *I, const int *J, const float *V, float *tempCSR, const float *vector_in, 
		float *vector_out, float *bp, float *pk, float *rk, int MAXIter)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for
	float dotp0,dotr0,dotr1,doth;

	//float *x=(float *) malloc(size1);
	int i,j;
	float threshold=0.0000001;
	int iter=0;
	float error,alphak,gamak, tempBp;
	int length,startCSR;
	error=1000;
	//initialize
	doth=0;
#pragma simd
	for (i=0;i<dimension;i++)
	{
		pk[i]=vector_in[i];
		rk[i]=vector_in[i];
		vector_out[i]=0;
		bp[i]=0;
		doth=doth+vector_in[i]*vector_in[i];
	}

	__assume_aligned(tempCSR,64);
	__assume_aligned(V,64);
	__assume_aligned(J,64);

	if (doth>0)
	{
		while (error>threshold&&iter<MAXIter){
			dotp0=0;
			dotr0=0;
			dotr1=0;
			//SpMV based on CSR format
			/*Prefetch the data: pk to L1 cache*/
			/*#pragma noprefetch
			  for (i=0;i<dimension;i+=16)
			  _mm_prefetch((const char*) &pk[i],0);*/

#pragma simd
			for (i=0;i<totalNum;i++)
				tempCSR[i]=V[i]*pk[J[i]];
			//#pragma noprefetch	
			for (i=0;i<dimension;i++)
			{	
				length=I[i+1]-I[i];
				startCSR=I[i];
				tempBp=0;
				/*Prefetch the V and J for next iterations*/
				/*for (j=I[i+1];j<I[i+2];j+=16)
				  _mm_prefetch((const char*) &V[j],0);*/
#pragma simd
				for(j=0;j<length;j++)
					tempBp+=tempCSR[j+startCSR];
				bp[i]=tempBp;
			}	
			//for (int k=0;k<5;k++) printf("pk %d at iter %d is %f\n",k, iter, pk[k]);
			//printf("CPU start\n");
			for (i=0;i<dimension;i++)
			{
				dotp0=dotp0+bp[i]*pk[i];
				dotr0=dotr0+rk[i]*rk[i];
			}	
			alphak=dotr0/dotp0;
#pragma simd
			for (i=0;i<dimension;i++) 
			{		
				vector_out[i]=vector_out[i]+alphak*pk[i];
				//if (i<10) printf ("pk[i] is %f, rk[i] is %f\n",pk[i],rk[i]);
				rk[i]=rk[i]-alphak*bp[i];
				dotr1=dotr1+rk[i]*rk[i];
			}
			gamak=dotr1/dotr0;
#pragma simd
			for (i=0;i<dimension;i++)
			{
				pk[i]=rk[i]+gamak*pk[i];
			}
			//printf("at iter %d, alphak is %f, gamak is %f\n",iter, alphak,gamak);
			error=sqrt(dotr1)/sqrt(doth);
			//printf("error at %d is %f\n",iter, error);
			iter++;
		}
	}
}

void fspai(S *SInput)
	/*void fspai(int *I, int *J, float *V, int *I_precond, int *J_precond, float *V_precond, const int maxRowNum, const int *numInRow, 
	  const int *rowNumAccum, const int *numInRowPrecond, const int *rowNumAccumPrecond, const float *diag, int colStart, int colEnd)*/
	/* totalNum: number of matrix, dimension: size of vector, I: row index, J: column index, V: matrix value, *_precond: preconditioner index and value
	*/
{
	int *I=SInput->I;
	int *J=SInput->J;
	float *V=SInput->V;
	int *I_precond=SInput->I_precond;
	int *J_precond=SInput->J_precond;
	float *V_precond=SInput->V_precond;
	int maxRowNum=SInput->maxRowNum;
	int *numInRow=SInput->numInRow;
	int *rowNumAccum=SInput->rowNumAccum;
	int *numInRowPrecond=SInput->numInRowPrecond;
	int *rowNumAccumPrecond=SInput->rowNumAccumPrecond;
	float *diag=SInput->diag;
	int colStart=SInput->colStart; 
	int colEnd=SInput->colEnd;
	int id=SInput->id;

	float *sortedV= (float *)malloc(maxRowNum*sizeof(float));
	int *sortedJ=(int *)malloc(maxRowNum*sizeof(int));

	int start,num;
	//printf("fspai start at %d end at %d maxRowNum is %d\n", colStart, colEnd,maxRowNum);
	//Sort the matrix, using insert sort algorithm, will try quick sort algorithm
	/*for (int i=colStart;i<colEnd;i++)
	  {

	  start=rowNumAccum[i];
	  num=numInRow[i];
	  insertSort(&J[start],&V[start],num,sortedJ,sortedV);
	  for (int j=0;j<num;j++)
	  {
	  J[start+j]=sortedJ[j];
	  V[start+j]=sortedV[j];
	  }

	  }*/
	//printf("fspai first finish\n");

	int subMatrixNum;
	float *subMatrix=(float *) _mm_malloc(maxRowNum*maxRowNum*sizeof(float),64);
	float *tempCSR=(float *) _mm_malloc(maxRowNum*maxRowNum*sizeof(float),64);
	int *subI=(int *) _mm_malloc((maxRowNum+1)*sizeof(int),64);
	int *subJ=(int *) _mm_malloc(maxRowNum*maxRowNum*sizeof(int),64);
	float *yk=(float *)malloc(maxRowNum*sizeof(float));
	float *xk=(float *)malloc(maxRowNum*sizeof(float));
	int *tempJ=(int *)malloc(maxRowNum*sizeof(int));
	float *bp=(float *) malloc(maxRowNum*sizeof(float));
	float *pk=(float *) malloc(maxRowNum*sizeof(float));
	float *rk=(float *) malloc(maxRowNum*sizeof(float));

	int iterNum;
	iterNum=100;
	float AY;
	float Lk;
	int subRowIndex, index1,index2,index3;
	int tempCol,colIndex;

	for (int i=colStart;i<colEnd;i+=1)
	{
		subRowIndex=0;
		start=0;
		subMatrixNum=0;
		for (int j=0;j<maxRowNum;j++)
		{
			xk[j]=0;
			yk[j]=0;
			subI[j]=0;
		}
		subI[maxRowNum]=0;
		/*Calculate the sub matrix and stored it as a CSR format*/
		for (int j=0; j<numInRow[i]; j++ )
		{
			index1=J[rowNumAccum[i]+j];
			if(index1>i)
			{
				//first time touch the column index larger than i
				if (start==0) start=j;
				colIndex=0;
				for (int k=0; k<numInRow[index1]; k++)
				{
					index2=rowNumAccum[index1]+k;
					if (J[index2]>i)
					{
						//travel at both row, the row i is currently in colIndex+j
						tempCol=J[rowNumAccum[i]+colIndex+start];
						//check the current line, whether it meet the 
						if (J[index2]==tempCol)
						{
							//subI[subMatrixNum]=subRowIndex;
							subI[subRowIndex+1]+=1;
							subJ[subMatrixNum]=colIndex;
							subMatrix[subMatrixNum]=V[index2];
							subMatrixNum++;
						}
						else 
						{
							while (J[index2]>tempCol&&colIndex+start<numInRow[i])
							{
								colIndex++;
								tempCol=J[rowNumAccum[i]+colIndex+start];
								if (J[index2]==tempCol)
								{
									//subI[subMatrixNum]=subRowIndex;
									subI[subRowIndex+1]+=1;
									subJ[subMatrixNum]=colIndex;
									subMatrix[subMatrixNum]=V[index2];
									subMatrixNum++;
								}
							}
						}
					}
				}
				xk[subRowIndex]=V[rowNumAccum[i]+j];
				tempJ[subRowIndex]=index1;
				subI[subRowIndex+1]+=subI[subRowIndex];
				subRowIndex=subRowIndex+1;
			}	
		}
		//calculate yk
		if (subMatrixNum>0)
		{
			solver(subRowIndex,subMatrixNum, subI, subJ, subMatrix, tempCSR, xk, yk, bp, pk, rk, iterNum);
		}
		//calculate Lk
		AY=0;
		//
		for (int p=0;p<subRowIndex;p++)
		{
			AY=AY+xk[p]*yk[p];
		}

		Lk = 1/sqrt(diag[i]-AY);
		//write the result
		for (int p=0;p<numInRowPrecond[i];p++)
		{
			index3=p+rowNumAccumPrecond[i];
			if (p==0)
			{
				I_precond[index3]=i;
				J_precond[index3]=i;
				V_precond[index3]=Lk;
			}
			else
			{
				I_precond[index3]=i;
				J_precond[index3]=tempJ[p-1];
				V_precond[index3]=(-Lk)*yk[p-1];
			}
		}
	}

}
