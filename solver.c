#include "kernel.h"
#include "solver.h"

void COO2ELL(const int *rowLocal, const int *colLocal, const float* matrixLocal, int **colELL,
	float **matrixELL, int **I_COO, int **J_COO, float **V_COO,const int *numInRow, 
	const int *rowNumAccum, const int localMatrixSize, const int localNumofRow, 
	int *sizeOut, int *max){

	int maxRowNum;
	maxRowNum=*max;
	int *hist=(int *)malloc((maxRowNum+1)*sizeof(int));
	int *accum=(int *)malloc((maxRowNum+1)*sizeof(int));

	int maxNew, sizeCOO=0, pointCOO=0;
	int rowBias=rowLocal[0];
	int numBias=rowNumAccum[rowBias];
	
	maxNew=ceil(((float) localMatrixSize)*1.2/((float) localNumofRow));
	for (int i=0;i<localNumofRow;i++)
	{
		int rowIndex=i+rowLocal[0];
		if (numInRow[i+rowBias]>maxNew)
			sizeCOO+=numInRow[i+rowBias]-maxNew;
	}	

	if (maxNew>0) maxRowNum=maxNew;
	
	*I_COO=(int *)malloc(sizeCOO*sizeof(int));
	*J_COO=(int *)malloc(sizeCOO*sizeof(int));
	*V_COO=(float *)malloc(sizeCOO*sizeof(float));
	*colELL=(int *)malloc(localNumofRow*maxRowNum*sizeof(int));
	*matrixELL=(float *)malloc(localNumofRow*maxRowNum*sizeof(float));
	
	int irregular=0;
	for (int i=0;i<localNumofRow;i++)
	{
		int rowIndex=i+rowLocal[0];
		//goto COO format
		if (numInRow[i+rowBias]>maxRowNum) {
			for (int j=0;j<numInRow[i+rowBias];j++){

				//the ELL value should still be set as zero
				if (j<maxRowNum){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[i]+j-numBias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
				}
				else{
					//assign COO value
					(*I_COO)[pointCOO]=rowLocal[rowNumAccum[rowIndex]+j-numBias];
					(*J_COO)[pointCOO]=colLocal[rowNumAccum[rowIndex]+j-numBias];
					(*V_COO)[pointCOO]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
					pointCOO++;
				}
			}
			irregular=irregular+1;
		}
		//goto ELL format
		else {
			for (int j=0;j<maxRowNum;j++){
				//write the ELL data
				if (j<numInRow[i+rowBias]){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[rowIndex]+j-numBias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
				}
				//write zero
				else{
					(*colELL)[i+j*localNumofRow]=0;
					(*matrixELL)[i+j*localNumofRow]=0;
				}
			}
		}
	}
	*max=maxRowNum;
	*sizeOut=sizeCOO;
	//printf("irregular row Number is %d\n",irregular);
	
}
void solver(const int dimension, const int totalNum, const int *I, const int *J, const float *V, const float *vector_in, 
			float *vector_out, float *error_track, int MAXIter)
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
	}
}

void solverGPU_HYB(const int dimension, const int totalNum, 
		const int *I_accum, const int* numInRow, 
		const int* I, const int *J, const float *V, 
		const int totalNumPrecond, 
		const int *I_precond_accum, const int *numInRowL, 
		const int *I_precond, const int *J_precond, const float *V_precond, 
		const int totalNumPrecondP,
		const int *I_precondP_accum, const int *numInRowLP, 
		const int *I_precondP, const int *J_precondP, const float *V_precondP, 
		const float *y, float *x,  
		const int MAXIter, int *realIter, const bool RODR, 
		const int partition_size, const int* part_boundary)
{

	int *colELL, *I_COO, *J_COO;
	float* matrixELL, *V_COO;
	int totalNumCOO, maxRowNum;
	int *col_d;
	float *V_d;
	int *I_COO_d, *J_COO_d;
	float *V_COO_d;
	COO2ELL(I,J,V,&colELL,&matrixELL,&I_COO, &J_COO, &V_COO, 
		numInRow, I_accum, totalNum, dimension, &totalNumCOO, &maxRowNum);
	cudaMalloc((void **) &col_d,maxRowNum*dimension*sizeof(int));
	cudaMalloc((void **) &V_d,maxRowNum*dimension*sizeof(float));
	cudaMemcpy(col_d,colELL,dimension*maxRowNum*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(V_d,matrixELL,dimension*maxRowNum*sizeof(float),cudaMemcpyHostToDevice);
	if (totalNumCOO>0){
		cudaMalloc((void **) &I_COO_d,totalNumCOO*sizeof(int));
		cudaMalloc((void **) &J_COO_d,totalNumCOO*sizeof(int));
		cudaMalloc((void **) &V_COO_d,totalNumCOO*sizeof(float));	
		cudaMemcpy(I_COO_d,I_COO,totalNumCOO*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(J_COO_d,J_COO,totalNumCOO*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(V_COO_d,V_COO,totalNumCOO*sizeof(float),cudaMemcpyHostToDevice);
	}	
/*matrix L'*/
	int *colELL_precond, *I_COO_L, *J_COO_L;
	float* matrixELL_precond, *V_COO_L;
	int totalNumCOO_L, maxRowNumPrecond;
	int *col_precond_d;
	float *V_precond_d;
	int *I_COO_L_d, *J_COO_L_d;
	float *V_COO_L_d;
	COO2ELL(I_precond,J_precond,V_precond,&colELL_precond, &matrixELL_precond, 
		&I_COO_L, &J_COO_L, &V_COO_L, numInRowL, I_precond_accum, totalNumPrecond, 
		dimension, &totalNumCOO_L, &maxRowNumPrecond);
	cudaMalloc((void **) &col_precond_d,maxRowNumPrecond*dimension*sizeof(int));
	cudaMalloc((void **) &V_precond_d,maxRowNumPrecond*dimension*sizeof(float));	
	cudaMemcpy(col_precond_d,colELL_precond,dimension*maxRowNumPrecond*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(V_precond_d,matrixELL_precond,dimension*maxRowNumPrecond*sizeof(float),cudaMemcpyHostToDevice);
	if (totalNumCOO_L>0){
		cudaMalloc((void **) &I_COO_L_d,totalNumCOO_L*sizeof(int));
		cudaMalloc((void **) &J_COO_L_d,totalNumCOO_L*sizeof(int));
		cudaMalloc((void **) &V_COO_L_d,totalNumCOO_L*sizeof(float));		
		cudaMemcpy(I_COO_L_d,I_COO_L,totalNumCOO_L*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(J_COO_L_d,J_COO_L,totalNumCOO_L*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(V_COO_L_d,V_COO_L,totalNumCOO_L*sizeof(float),cudaMemcpyHostToDevice);
	}	
	/*matrix L*/
	int *colELL_precondP, *I_COO_LP, *J_COO_LP;
	float* matrixELL_precondP, *V_COO_LP;
	int totalNumCOO_LP, maxRowNumPrecondP;
	int *col_precondP_d;
	float *V_precondP_d;
	int *I_COO_LP_d, *J_COO_LP_d;
	float *V_COO_LP_d;

	COO2ELL(I_precondP,J_precondP,V_precondP,&colELL_precondP, &matrixELL_precondP, 
		&I_COO_LP, &J_COO_LP, &V_COO_LP, numInRowLP, I_precondP_accum, totalNumPrecond, 
		dimension, &totalNumCOO_LP, &maxRowNumPrecondP);

	cudaMalloc((void **) &col_precondP_d,maxRowNumPrecondP*dimension*sizeof(int));
	cudaMalloc((void **) &V_precondP_d,maxRowNumPrecondP*dimension*sizeof(float));	
	cudaMemcpy(col_precondP_d,colELL_precondP,dimension*maxRowNumPrecondP*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(V_precondP_d,matrixELL_precondP,dimension*maxRowNumPrecondP*sizeof(float),cudaMemcpyHostToDevice);
	if (totalNumCOO_LP>0){
		cudaMalloc((void **) &I_COO_LP_d,totalNumCOO_LP*sizeof(int));
		cudaMalloc((void **) &J_COO_LP_d,totalNumCOO_LP*sizeof(int));
		cudaMalloc((void **) &V_COO_LP_d,totalNumCOO_LP*sizeof(float));		
		cudaMemcpy(I_COO_LP_d,I_COO_LP,totalNumCOO_LP*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(J_COO_LP_d,J_COO_LP,totalNumCOO_LP*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(V_COO_LP_d,V_COO_LP,totalNumCOO_LP*sizeof(float),cudaMemcpyHostToDevice);
	}

	size_t size0=dimension*sizeof(float);
	float *rk_d;//r0 and r1
	float *pk_d;//p0 and p1
	float *bp_d;
	float *zk_d;
	float *zk1_d;
	float *x_d;
	int* part_boundary_d;
	
	float *rk=(float *)malloc(size0);
	float *zk=(float *)malloc(size0);
	float *zk_1=(float *)malloc(size0);
	
	float *vector_in_d;
	float alphak,betak,dotp0, dotrz0,dotrz1,doth,alphak_1,dotr0;
	
	//float dotz0,dotz0_compare;
	
	float error=10000;
	float threshold;
	cudaMalloc((void **) &rk_d,size0);
	cudaMalloc((void **) &pk_d,size0);
	cudaMalloc((void **) &zk_d,size0);
	cudaMalloc((void **) &zk1_d,size0);
	cudaMalloc((void **) &bp_d,size0);
	cudaMalloc((void **) &x_d,size0);
	cudaMalloc((void **) &part_boundary_d, partition_size*sizeof(int));

	cudaMalloc((void **) &vector_in_d,size0);
	//cudaMalloc((void **) &vector_out_d,size0);
	
	cudaMemcpy(vector_in_d, vector_in, size0, cudaMemcpyHostToDevice);	
	cudaMemcpy(part_boundary_d, part_boundaryd, partition_size*sizeof(int), cudaMemcpyHostToDevice);	
	
	//
	initialize_bp(dimension,zk_d);
	initialize_bp(dimension,zk1_d);
	initialize_r(dimension, rk_d, vector_in_d);
		

	
	cublasHandle_t handle;
	cublasCreate(&handle);			
	if(RODR){
		matrix_vectorELL(dimension, dimension, maxRowNumPrecond, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
				true, part_boundary_d);
	}
	else{
		matrix_vectorELL(dimension, dimension, maxRowNumPrecond, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
				false, NULL);
	}
	if (totalNumCOO_L>0) matrix_vectorCOO(totalNumCOO_L, I_COO_L, J_COO_L, V_COO_L, rk_d, zk1_d, 0, 0);
	if(RODR){
		matrix_vectorELL(dimension, dimension, maxRowNumPrecondP, col_precondP_d, 
			V_precondP_d, zk1_d, zk_d, 0, 0, true, part_boundary_d);
	}
	else{
		matrix_vectorELL(dimension, dimension, maxRowNumPrecondP, col_precondP_d, 
			V_precondP_d, zk1_d, zk_d, 0, 0, false, NULL);
	}
	
	if (totalNumCOO_LP>0) matrix_vectorCOO(totalNumCOO_LP, I_COO_LP, J_COO_LP, V_COO_LP, zk1_d, zk_d,0,0);
	
	//cublasSdot(handle,dimension,zk_d,1,zk_d,1,&dotz0);
	//for (int i=0;i<5;i++) printf("dotz0 is %f, dotz0_compare is %f\n",dotz0,dotz0_compare);
	//give an initial value as r0, p0, and set all 
	initialize_all(dimension,pk_d,bp_d,x_d,zk_d,vector_in_d);
	
	
	cublasSdot(handle, dimension,vector_in_d,1,vector_in_d,1,&doth);
	float errorNorm=sqrt(doth);
	int iter=0;
	
	while (iter<MAXIter&&error>1E-7){
		//matrix-vector multiplication, accelerated by our own functions
		if(RODR){
			matrix_vectorELL(dimension, dimension, maxRowNum, col_d,V_d,pk_d,bp_d,0,0,
			true, part_boundary_d);
		}
		else{
			matrix_vectorELL(dimension, dimension, maxRowNum, col_d,V_d,pk_d,bp_d,0,0,
			false, NULL);
		}
		if (totalNumCOO>0) matrix_vectorCOO(totalNumCOO, I_COO, J_COO, V_COO, pk_d, bp_d, 0, 0);
		//vector-vectror multiplication
		cublasSdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		//r_k*z_k
		cublasSdot(handle,dimension,rk_d,1,zk_d,1,&dotrz0);
		/*cudaMemcpy(rk, rk_d, size0, cudaMemcpyDeviceToHost);
		cudaMemcpy(zk, zk_d, size0, cudaMemcpyDeviceToHost);
		float dotrz0_compare=0;
		for (int k=0;k<dimension;k++)
		{
			dotrz0_compare=dotrz0_compare+rk[k]*zk[k];
		}		
		printf("dotrz0 is %f, dotrz0_compare is %f\n",dotrz0,dotrz0_compare);*/
		alphak=dotrz0/dotp0;
		alphak_1=-alphak;
		
		//update x_k to x_k+1, x=x+alphak*p;
		cublasSaxpy(handle,dimension,&alphak,pk_d,1,x_d,1);
		//update r_k to r_k+1, r=r-alphak*p;
		cublasSaxpy(handle,dimension,&alphak_1,bp_d,1,rk_d,1);
		initialize_bp(dimension,zk_d);
		initialize_bp(dimension,zk1_d);
		if(RODR){
			matrix_vectorELL(dimension, dimension, maxRowNumPrecond, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
			true, part_boundary_d);
		}
		else{
			matrix_vectorELL(dimension, dimension, maxRowNumPrecond, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
			false, NULL);
		}
		if (totalNumCOO_L>0) matrix_vectorCOO(totalNumCOO_L, I_COO_L, J_COO_L, V_COO_L, rk_d,zk1_d,0 ,0);
		if(RODR){	
			matrix_vectorell(dimension, dimension, maxrownumprecondp, col_precondp_d,v_precondp_d,
				zk1_d,zk_d,0,0, true, part_boundary_d);
		}
		else{
			matrix_vectorell(dimension, dimension, maxrownumprecondp, col_precondp_d,v_precondp_d,
				zk1_d,zk_d,0,0, false, NULL);
		}
		if (totalNumCOO_LP>0) matrix_vectorCOO(totalNumCOO_LP, I_COO_LP, J_COO_LP, V_COO_LP, zk1_d,zk_d,0,0 );
		//r_k+1 * z_k+1
		cublasSdot(handle,dimension,rk_d,1,zk_d,1,&dotrz1);
		betak=dotrz1/dotrz0;
		//printf("dotp0 is %f, dotrz0 is %f dotrz1 is %f, betak is %f and alphak is %f at iter %d\n", dotp0, dotrz0,dotrz1,betak, alphak, iter);
		//p=r+gamak*p;
		myxpy(dimension,betak,zk_d,pk_d);//if using cublas, we need two functin to complete it
		cublasSdot(handle,dimension,rk_d,1,rk_d,1,&dotr0);
		error=sqrt(dotr0)/errorNorm;
		initialize_bp(dimension,bp_d);
		cudaDeviceSynchronize();
		error_track[iter]=error;
		
		//printf("at iter %d, alpha is %f, beta is %f, dotrz0 is %f, dotrz1 is %f, dotp0 is %f\n", iter, alphak, betak, dotrz0, dotrz1, dotp0);
		iter++;
	}
	printf("max iteration is %d with error %e\n",iter, error);
	cublasSestroy(handle);
	cudaMemcpy(vector_out, x_d, size0, cudaMemcpyDeviceToHost);
	free(I_COO);
	free(J_COO);
	free(V_COO);
	free(matrixELL);
	free(colELL);
}

void solverPrecondCPU(const int dimension, const int totalNum, const int *I_accum, const int *J, 
		const float *V, const int totalNumPrecond, const int *I_precond_accum, 
		const int *J_precond, const float *V_precond, const int totalNumPrecondP,
		const int *I_precondP_accum, const int *J_precondP, const float *V_precondP, 
		const float *vector_in, float *vector_out, const int MAXIter, int *realIter)
{
	
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

	stride=ceil((float)totalNum/procNum);
	stridePrecond=ceil((float)totalNumPrecond/procNum);
	stridePrecondP=ceil((float)totalNumPrecondP/procNum);
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


	/*#pragma omp parallel for
	  for (int i=0;i<totalNumPrecond;i++)
	  zk1[I_precond[i]]+=V_precond[i]*rk[J_precond[i]];*/
#pragma omp parallel private(rank,tempV)
	{
		rank=omp_get_thread_num();
		for (int i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
			tempV=0;
			a=I_precond_accum[i];
			b=I_precond_accum[i+1];
#pragma simd
			for (int j=a;j<b;j++){
				tempV+=V_precond[j]*rk[J_precond[j]];
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
					tempV+=V[j]*pk[J[j]];
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
					tempV+=V_precond[j]*rk[J_precond[j]];
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
	/*--------------finish the iterations---------------------*/
	gettimeofday(&end1, NULL);	
	printf("max iteration is %d with error %e\n",iter, error);

	/*------------------------------------------------------------Finish the iterative solver-------------------------------------------------*/



	gettimeofday(&end1, NULL);	
	printf("Solver Xeon_phi time is %ld us\n",(end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec));
	float timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, Xeon_Phi Gflops is %f\n ",iter, (1e-9*(totalNum*4+14*dimension)*1000*iter)/timeByMs);



}


