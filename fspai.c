#include "fspai.h"

void formatChange(int dimension, int *numInRow, int *totalNum_1, int *I, int *J, double *V, int **rowNumAccum_1, int **J_1, double **V_1)
{
	int totalNum;
	int *numInRow_1=(int *)malloc((dimension)*sizeof(int));
	*rowNumAccum_1=(int *)malloc((dimension+1)*sizeof(int));
	int *rowNumAccum=(int *)malloc((dimension+1)*sizeof(int));
	(*rowNumAccum_1)[0]=0;
	rowNumAccum[0]=0;
	for (int i=0;i<dimension;i++){
		numInRow_1[i]=ceil((double)numInRow[i]/padding)*padding;
		rowNumAccum[i+1]=rowNumAccum[i]+numInRow[i];
		(*rowNumAccum_1)[i+1]=(*rowNumAccum_1)[i]+numInRow_1[i];
	}
	//for (int i=0;i<10;i++) printf("%d\t%d\n",numInRow[i],numInRow_1[i]);
	
	totalNum=(*rowNumAccum_1)[dimension];
	
	*J_1=(int *)_mm_malloc(totalNum*sizeof(int),64);
	*V_1=(double *)_mm_malloc(totalNum*sizeof(double),64);
	//printf("totalNum is %d\n",totalNum);
	
	#pragma omp parallel for
	for (int i=0;i<dimension;i++){
		for (int j=0;j<numInRow_1[i];j++){
			
			if (j<numInRow[i]){
				(*J_1)[(*rowNumAccum_1)[i]+j]=J[rowNumAccum[i]+j];
				(*V_1)[(*rowNumAccum_1)[i]+j]=V[rowNumAccum[i]+j];
			}
			else
			{
				(*J_1)[(*rowNumAccum_1)[i]+j]=J[rowNumAccum[i]+numInRow[i]-1];
				(*V_1)[(*rowNumAccum_1)[i]+j]=0;
			}
		}
	}
	*totalNum_1=totalNum;
}
/*double for preconditioner generating only*/
void solverCPU(const int dimension, const int totalNum, const int *I, const int *J, const float*V, const float *vector_in, 
			float *vector_out, float *error_track, int MAXIter, int *realIter)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for
	float dotp0,dotr0,dotr1,doth;
	size_t size1=dimension*sizeof(float);
	float *bp=(float *) malloc(size1);
	float *pk=(float *) malloc(size1);
	float *rk=(float *) malloc(size1);
	//float *x=(float *) malloc(size1);
	int i;
	float threshold=0.0000001;
	int iter=0;
	float error,alphak,gamak;
	error=1000;
	//initialize
	doth=0;
	for (i=0;i<dimension;i++)
	{
		pk[i]=vector_in[i];
		rk[i]=vector_in[i];
		vector_out[i]=0;
		bp[i]=0;
		doth=doth+vector_in[i]*vector_in[i];
	}
	*realIter=0;
	if (doth>0)
	{
		while (error>threshold&&iter<MAXIter){
			dotp0=0;
			dotr0=0;
			dotr1=0;
			for (i=0;i<totalNum;i++)
			{		
				bp[I[i]]+=V[i]*pk[J[i]];//multiplication
			}	
			//for (int k=0;k<5;k++) printf("pk %d at iter %d is %f\n",k, iter, pk[k]);
			//printf("CPU start\n");
			for (i=0;i<dimension;i++)
			{
				dotp0=dotp0+bp[i]*pk[i];
				dotr0=dotr0+rk[i]*rk[i];
			}	
			alphak=dotr0/dotp0;
			for (i=0;i<dimension;i++) 
			{		
				vector_out[i]=vector_out[i]+alphak*pk[i];
				//if (i<10) printf ("pk[i] is %f, rk[i] is %f\n",pk[i],rk[i]);
				rk[i]=rk[i]-alphak*bp[i];
				dotr1=dotr1+rk[i]*rk[i];
			}
			gamak=dotr1/dotr0;
			for (i=0;i<dimension;i++)
			{
				pk[i]=rk[i]+gamak*pk[i];
				bp[i]=0;
			}
			//printf("at iter %d, alphak is %f, gamak is %f\n",iter, alphak,gamak);
			error=sqrt(dotr1)/sqrt(doth);
			error_track[iter]=error;
			//printf("error at %d is %f\n",iter, error);
			iter++;
		}
		*realIter=iter;
	}
}

void solverPrecondCOO(const int dimension, const int totalNum, const int *I_accum, const int *J, const double *V, const int totalNumPrecond, const int *I_precond_accum,
				const int *J_precond, const double *V_precond, const int totalNumPrecondP, const int *I_precondP_accum, const int *J_precondP, const double *V_precondP,
				const double *y, double *x, const int MAXIter, int *realIter)
{
	size_t size0=dimension*sizeof(double);
	
	double *rk=(double *)malloc(size0);
	double *zk=(double *)malloc(size0);
	double *zk1=(double *)malloc(size0);
	double *pk=(double *)malloc(size0);
	double *bp=(double *)malloc(size0);
	
	double *vector_in_d;
	double alphak,betak,dotp0, dotrz0,dotrz1,doth,alphak_1,dotr0;
	double alpha = 1.0f, beta = 0.0f; 
	
	//double dotz0,dotz0_compare;
	
	double error=10000;
	double threshold;

	for (int i=0;i<dimension;i++){
		zk[i]=0;
		zk1[i]=0;
		rk[i]=y[i];
	}
	
	int procNum, rank;
	double tempV;
	
	#pragma omp parallel
	{
		#pragma omp master
		{
			procNum=omp_get_num_threads();
		}
	}
	//printf("thread num is %d\n",procNum);
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
	
	stride=totalNum/procNum;
	stridePrecond=totalNumPrecond/procNum;
	stridePrecondP=totalNumPrecondP/procNum;
	int bias, biasPrecond,biasPrecondP;
	
	bias=1;
	biasPrecond=1;
	biasPrecondP=1;

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
	//for (int i=0;i<=6;i++) printf("at i=%d, bias is %d, biasPrecond is %d, biasPrecondP is %d\n",i, boundary[i],precond_boundary[i],precondP_boundary[i]);
	char *matdescra="GNNC  ";
	/*#pragma omp parallel for
	for (int i=0;i<totalNumPrecond;i++)
			zk1[I_precond[i]]+=V_precond[i]*rk[J_precond[i]];*/
	#pragma omp parallel private(rank,tempV)
	{
		rank=omp_get_thread_num();
		for (int i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
			tempV=0;
			#pragma simd
			for (int j=I_precond_accum[i];j<I_precond_accum[i+1];j++){
				tempV+=V_precond[j]*rk[J_precond[j]];
				//#pragma omp master
				//if (i<8) printf("V is %f, J is %d, rk is %f, tempv is %f\n",V_precond[i+j], J_precond[i+j], rk[J_precond[i+j]], tempV);
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
			#pragma simd
			for (int j=I_precondP_accum[i];j<I_precondP_accum[i+1];j++)
				tempV+=V_precondP[j]*zk1[J_precondP[j]];
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
	double errorNorm=sqrt(doth);
	printf("errorNorm is %f\n",errorNorm);
	int iter=0;
	
	while (iter<MAXIter&&error>1E-7){
		//matrix-vector multiplication, accelerated by our own functions
		/*#pragma omp parallel for
		for (int i=0;i<totalNum;i++)
			bp[I[i]]+=V[i]*pk[J[i]];*/
			
		#pragma omp parallel private(rank,tempV)
		{
			rank=omp_get_thread_num();
			for (int i=boundary[rank];i<boundary[rank+1];i++){
				tempV=0;
				#pragma simd
				for (int j=I_accum[i];j<I_accum[i+1];j++){
					tempV+=V[j]*pk[J[j]];
					//#pragma omp master
					//if (i<8) printf("V is %f, J is %d, rk is %f, tempv is %f\n",V_precond[i+j], J_precond[i+j], rk[J_precond[i+j]], tempV);
				}
				bp[i]+=tempV;
				//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
			}
		}
		//vector-vectror multiplication
		dotp0=0;
		#pragma omp parallel for reduction(+:dotp0)
		for (int i=0;i<dimension;i++)
			dotp0+=bp[i]*pk[i];
		
		//r_k*z_k
		//cublasSdot(handle,dimension,rk_d,1,zk_d,1,&dotrz0);
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
		
		//mkl_scoomv("N", &dimension, &dimension, &alpha, matdescra, V_precond, I_precond, J_precond, &totalNumPrecond, rk, &beta, zk1);
		/*for (int i=0;i<totalNum;i+=8){
			#pragma simd
			for (int j=0;j<8;j++)
				tempV+=;
			zk1[]=tempV;
		}*/
		#pragma omp parallel private(rank,tempV)
		{
			rank=omp_get_thread_num();
			for (int i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
				tempV=0;
				#pragma simd
				for (int j=I_precond_accum[i];j<I_precond_accum[i+1];j++){
					tempV+=V_precond[j]*rk[J_precond[j]];
					//#pragma omp master
					//if (i<8) printf("V is %f, J is %d, rk is %f, tempv is %f\n",V_precond[i+j], J_precond[i+j], rk[J_precond[i+j]], tempV);
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
				#pragma simd
				for (int j=I_precondP_accum[i];j<I_precondP_accum[i+1];j++)
					tempV+=V_precondP[j]*zk1[J_precondP[j]];
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
	printf("max iteration is %d with error %e\n",iter, error);
	*realIter=iter;
}


void solverPrecondPhi(const int dimension, const int totalNum, const int *I, const int *J, const double *V, const int totalNumPrecond, const int *I_precond,
				const int *J_precond, const double *V_precond, const int totalNumPrecondP, const int *I_precondP, const int *J_precondP, const double *V_precondP,
				const double *y, double *x, const int MAXIter, int *realIter, int rank)
{
	size_t size0=dimension*sizeof(double);
	
	double *rk=(double *)malloc(size0);
	double *zk=(double *)malloc(size0);
	double *zk1=(double *)malloc(size0);
	double *pk=(double *)malloc(size0);
	double *bp=(double *)malloc(size0);
	
	double *vector_in_d;
	double alphak,betak,dotp0, dotrz0,dotrz1,doth,alphak_1,dotr0;
	double alpha = 1.0f, beta = 0.0f; 
	
	//double dotz0,dotz0_compare;
	
	double error=10000;
	double threshold;

	for (int i=0;i<dimension;i++){
		zk[i]=0;
		zk1[i]=0;
		rk[i]=y[i];
	}
	
	int procNum;
	double tempV;
	
	#pragma omp parallel
	{
		#pragma omp master
		{
			procNum=omp_get_num_threads();
		}
	}
	//printf("thread num is %d\n",procNum);
	int *boundary=(int *)malloc(procNum+1*sizeof(int));
	int *precond_boundary=(int *)malloc(procNum+1*sizeof(int));
	int *precondP_boundary=(int *)malloc(procNum+1*sizeof(int));
	boundary[0]=0;
	precond_boundary[0]=0;
	precondP_boundary[0]=0;
	boundary[procNum]=totalNum;
	precond_boundary[procNum]=totalNumPrecond;
	precondP_boundary[procNum]=totalNumPrecondP;
	
	int bias;
	for (int i=1;i<procNum;i++){
		bias=(totalNum/procNum)*i;
		while(I[bias]==I[bias+1]) bias++;
		boundary[i]=bias+1;
		
		/*Determine the boundary of precond matrix, reset the value of bias*/
		bias=(totalNumPrecond/procNum)*i;
		while(I_precond[bias]==I_precond[bias+1]) bias++;
		precond_boundary[i]=bias+1;
		
		/*Determine the boundary of precondP matrix, reset the value of bias*/
		bias=(totalNumPrecondP/procNum)*i;
		while(I_precondP[bias]==I_precondP[bias+1]) bias++;
		precondP_boundary[i]=bias+1;
		
		//printf("in %d, boundary is %d, precond_boundary is %d, precondP_boundary is %d\n",i,boundary[i],precond_boundary[i],precondP_boundary[i] );
		
	}
	
	char *matdescra="GNNC  ";
	/*#pragma omp parallel for
	for (int i=0;i<totalNumPrecond;i++)
			zk1[I_precond[i]]+=V_precond[i]*rk[J_precond[i]];*/
	#pragma omp parallel private(rank,tempV)
	{
		rank=omp_get_thread_num();
		for (int i=precond_boundary[rank];i<precond_boundary[rank+1];i+=padding){
			tempV=0;
			#pragma simd
			for (int j=0;j<padding;j++){
				tempV+=V_precond[i+j]*rk[J_precond[i+j]];
				//#pragma omp master
				//if (i<8) printf("V is %f, J is %d, rk is %f, tempv is %f\n",V_precond[i+j], J_precond[i+j], rk[J_precond[i+j]], tempV);
			}
			zk1[I_precond[i]]+=tempV;
			//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
		}
	}
	
	rank=0;
	#pragma omp parallel private(rank,tempV)
	{
		rank=omp_get_thread_num();
		for (int i=precondP_boundary[rank];i<precondP_boundary[rank+1];i++)
				zk[I_precondP[i]]+=V_precondP[i]*zk1[J_precondP[i]];
	}
			
			
	/*#pragma omp parallel private(rank,tempV)
	{
		rank=omp_get_thread_num();
		for (int i=precondP_boundary[rank];i<precondP_boundary[rank+1];i+=8){
			tempV=0;
			#pragma simd
			for (int j=0;j<8;j++)
				tempV+=V_precondP[i+j]*zk1[J_precondP[i+j]];
			zk[I_precondP[i]]+=tempV;
		}
		
	}*/
	
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
	double errorNorm=sqrt(doth);
	//printf("errorNorm is %f\n",errorNorm);
	int iter=0;
	
	while (iter<MAXIter&&error>1E-7){
		//matrix-vector multiplication, accelerated by our own functions
		//mkl_scoomv("N", &dimension, &dimension, &alpha, matdescra, V, I, J, &totalNum, pk, &beta, bp);
		#pragma omp parallel for
		for (int i=0;i<totalNum;i++)
			bp[I[i]]+=V[i]*pk[J[i]];
		//vector-vectror multiplication
		//cublasSdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		dotp0=0;
		#pragma omp parallel for reduction(+:dotp0)
		for (int i=0;i<dimension;i++)
			dotp0+=bp[i]*pk[i];
		
		//r_k*z_k
		//cublasSdot(handle,dimension,rk_d,1,zk_d,1,&dotrz0);
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
		
		//mkl_scoomv("N", &dimension, &dimension, &alpha, matdescra, V_precond, I_precond, J_precond, &totalNumPrecond, rk, &beta, zk1);
		/*for (int i=0;i<totalNum;i+=8){
			#pragma simd
			for (int j=0;j<8;j++)
				tempV+=;
			zk1[]=tempV;
		}*/
		#pragma omp parallel for
		for (int i=0;i<totalNumPrecond;i++)
			zk1[I_precond[i]]+=V_precond[i]*rk[J_precond[i]];
		//mkl_scoomv("N", &dimension, &dimension, &alpha, matdescra, V_precondP, I_precondP, J_precondP, &totalNumPrecondP, zk1, &beta, zk);
		#pragma omp parallel for
		for (int i=0;i<totalNumPrecondP;i++){			
			zk[I_precondP[i]]+=V_precondP[i]*zk1[J_precondP[i]];
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
	printf("max iteration is %d with error %e\n",iter, error);
	*realIter=iter;
}


void fspaiCPU(S *SInput)
/*void fspai(int *I, int *J, double *V, int *I_precond, int *J_precond, double *V_precond, const int maxRowNum, const int *numInRow, 
			const int *rowNumAccum, const int *numInRowPrecond, const int *rowNumAccumPrecond, const double *diag, int colStart, int colEnd)*/
/* totalNum: number of matrix, dimension: size of vector, I: row index, J: column index, V: matrix value, *_precond: preconditioner index and value
*/
{
	int *I=SInput->I;
	int *J=SInput->J;
	double *V=SInput->V;
	int *I_precond=SInput->I_precond;
	int *J_precond=SInput->J_precond;
	double *V_precond=SInput->V_precond;
	int maxRowNum=SInput->maxRowNum;
	int *numInRow=SInput->numInRow;
	int *rowNumAccum=SInput->rowNumAccum;
	int *numInRowPrecond=SInput->numInRowPrecond;
	int *rowNumAccumPrecond=SInput->rowNumAccumPrecond;
	double *diag=SInput->diag;
	int colStart=SInput->colStart; 
	int colEnd=SInput->colEnd;
	int id=SInput->id;
	
	double *sortedV= (double *)malloc(maxRowNum*sizeof(double));
	int *sortedJ=(int *)malloc(maxRowNum*sizeof(int));
	
	int start,num;
	//printf("fspai start at %d end at %d maxRowNum is %d\n", colStart, colEnd,maxRowNum);
	//Sort the matrix, using insert sort algorithm, will try quick sort algorithm
	for (int i=colStart;i<colEnd;i++)
	{
		
		start=rowNumAccum[i];
		num=numInRow[i];
		insertSort(&J[start],&V[start],num,sortedJ,sortedV);
		for (int j=0;j<num;j++)
		{
			J[start+j]=sortedJ[j];
			V[start+j]=sortedV[j];
		}
		
	}
	//printf("fspai first finish\n");
	
	int subMatrixNum;
	float *subMatrix=(float *) malloc(maxRowNum*maxRowNum*sizeof(double));
	int *subI=(int *) malloc(maxRowNum*maxRowNum*sizeof(int));
	int *subJ=(int *) malloc(maxRowNum*maxRowNum*sizeof(int));
	float *yk=(float *)malloc(maxRowNum*sizeof(float));
	float *xk=(float *)malloc(maxRowNum*sizeof(float));
	int *tempJ=(int *)malloc(maxRowNum*sizeof(int));
	
	int iterNum;
	iterNum=100;
	float *error=(float *)malloc(iterNum*sizeof(float));
	float AY;
	float Lk;
	int subRowIndex, index1,index2,index3;
	int tempCol,colIndex;
	//Complete computation in each column seperately
	for (int i=colStart;i<colEnd;i++)
	{
		subRowIndex=0;
		start=0;
		subMatrixNum=0;
		for (int j=0;j<maxRowNum;j++)
		{
			xk[j]=0;
			yk[j]=0;
		}
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
							subI[subMatrixNum]=subRowIndex;
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
									subI[subMatrixNum]=subRowIndex;
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
				subRowIndex=subRowIndex+1;
			}	
		}
		
		//calculate yk
		int realIter;
		if (subMatrixNum>0)
		{
			solverCPU(subRowIndex,subMatrixNum, subI, subJ, subMatrix, xk, yk, error, iterNum, &realIter);
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
				if (abs(V_precond[index3])>10000)
					printf("error happend at line %d colum %d of preconditioner\n", i, tempJ[p-1]);
			}
		}
	}

}
