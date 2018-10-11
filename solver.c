#include "kernel.h"
#include "solver.h"
#include "convert.h"

void solver(const unsigned int dimension, const unsigned int totalNum, const unsigned int *I, const unsigned int *J, const double *V, const double *vector_in, 
			double *vector_out, double *error_track, unsigned int MAXIter)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for
	double dotp0,dotr0,dotr1,doth;
	size_t size1=dimension*sizeof(double);
	double *bp=(double *) malloc(size1);
	double *pk=(double *) malloc(size1);
	double *rk=(double *) malloc(size1);
	//double *x=(double *) malloc(size1);
	unsigned int i;
	double threshold=0.0000001;
	unsigned int iter=0;
	double error,alphak,gamak;
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
	if (doth>0){
		while (error>threshold&&iter<MAXIter){
			dotp0=0;
			dotr0=0;
			dotr1=0;
			for (i=0;i<totalNum;i++)
			{		
				bp[I[i]]+=V[i]*pk[J[i]];//multiplication
			}	
			//for (unsigned int k=0;k<5;k++) printf("pk %d at iter %d is %f\n",k, iter, pk[k]);
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

void solverGPU_HYB(const unsigned int dimension, 
		const unsigned int totalNum, const unsigned int* numInRow, unsigned int maxRowNum,
		const unsigned int *row_idx,  const unsigned int* I, const unsigned int *J, const double *V, 
		const unsigned int totalNumPrecond, const unsigned int *numInRowL, unsigned int maxRowNumPrecond,
		const unsigned int *row_idxL,  
		const unsigned int *I_precond, const unsigned int *J_precond, const double *V_precond, 
		const unsigned int totalNumPrecondP, const unsigned int *numInRowLP, unsigned int maxRowNumPrecondP,
		const unsigned int *row_idxLP,  
		const unsigned int *I_precondP, const unsigned int *J_precondP, const double *V_precondP, 
		const double *vector_in, double *vector_out,  
		const unsigned int MAXIter, unsigned int *realIter,  const cb_s cb,
		const unsigned int partition_size, const unsigned int* part_boundary)
{

	unsigned int *colELL, *I_COO, *J_COO, ELL_width, ELL_widthL, ELL_widthLP;
	double* matrixELL, *V_COO;
	unsigned int totalNumCOO; 
	unsigned int *col_d;
	double *V_d;
	bool RODR, BOLCK;
	unsigned int *I_COO_d, *J_COO_d;
	double *V_COO_d;
	RODR = cb.RODR;
	BLOCK = cb.BLOCK;
	if(BLOCK){
		
		COO2ELL_block();
	} else {
		COO2ELL(I,J,V,&colELL,&matrixELL,&I_COO, &J_COO, &V_COO, 
			numInRow, row_idx, totalNum, dimension, &totalNumCOO, maxRowNum, &ELL_width);
	}
	
	cudaMalloc((void **) &col_d,ELL_width*dimension*sizeof(unsigned int));
	cudaMalloc((void **) &V_d,ELL_width*dimension*sizeof(double));
	cudaMemcpy(col_d,colELL,dimension*ELL_width*sizeof(unsigned int),cudaMemcpyHostToDevice);
	cudaMemcpy(V_d,matrixELL,dimension*ELL_width*sizeof(double),cudaMemcpyHostToDevice);
	if(BLOCK){
		
	}
	printf("ELL_width is %d, and totalNumCOO is %d\n", ELL_width, totalNumCOO);
	if (totalNumCOO>0){
		cudaMalloc((void **) &I_COO_d,totalNumCOO*sizeof(unsigned int));
		cudaMalloc((void **) &J_COO_d,totalNumCOO*sizeof(unsigned int));
		cudaMalloc((void **) &V_COO_d,totalNumCOO*sizeof(double));	
		cudaMemcpy(I_COO_d,I_COO,totalNumCOO*sizeof(unsigned int),cudaMemcpyHostToDevice);
		cudaMemcpy(J_COO_d,J_COO,totalNumCOO*sizeof(unsigned int),cudaMemcpyHostToDevice);
		cudaMemcpy(V_COO_d,V_COO,totalNumCOO*sizeof(double),cudaMemcpyHostToDevice);
	}	
/*matrix L'*/
	unsigned int *colELL_precond, *I_COO_L, *J_COO_L;
	double* matrixELL_precond, *V_COO_L;
	unsigned int totalNumCOO_L;
	unsigned int *col_precond_d;
	double *V_precond_d;
	unsigned int *I_COO_L_d, *J_COO_L_d;
	double *V_COO_L_d;
	COO2ELL(I_precond,J_precond,V_precond,&colELL_precond, &matrixELL_precond, 
		&I_COO_L, &J_COO_L, &V_COO_L, numInRowL, row_idxL, totalNumPrecond, 
		dimension, &totalNumCOO_L, maxRowNumPrecond, &ELL_widthL);
	cudaMalloc((void **) &col_precond_d,ELL_widthL*dimension*sizeof(unsigned int));
	cudaMalloc((void **) &V_precond_d,ELL_widthL*dimension*sizeof(double));	
	cudaMemcpy(col_precond_d,colELL_precond,dimension*ELL_widthL*sizeof(unsigned int),cudaMemcpyHostToDevice);
	cudaMemcpy(V_precond_d,matrixELL_precond,dimension*ELL_widthL*sizeof(double),cudaMemcpyHostToDevice);
	printf("ELL_widthL is %d, and totalNumCOO_L is %d\n", ELL_widthL, totalNumCOO_L);
	if (totalNumCOO_L>0){
		cudaMalloc((void **) &I_COO_L_d,totalNumCOO_L*sizeof(unsigned int));
		cudaMalloc((void **) &J_COO_L_d,totalNumCOO_L*sizeof(unsigned int));
		cudaMalloc((void **) &V_COO_L_d,totalNumCOO_L*sizeof(double));		
		cudaMemcpy(I_COO_L_d,I_COO_L,totalNumCOO_L*sizeof(unsigned int),cudaMemcpyHostToDevice);
		cudaMemcpy(J_COO_L_d,J_COO_L,totalNumCOO_L*sizeof(unsigned int),cudaMemcpyHostToDevice);
		cudaMemcpy(V_COO_L_d,V_COO_L,totalNumCOO_L*sizeof(double),cudaMemcpyHostToDevice);
	}	
	/*matrix L*/
	unsigned int *colELL_precondP, *I_COO_LP, *J_COO_LP;
	double* matrixELL_precondP, *V_COO_LP;
	unsigned int totalNumCOO_LP;
	unsigned int *col_precondP_d;
	double *V_precondP_d;
	unsigned int *I_COO_LP_d, *J_COO_LP_d;
	double *V_COO_LP_d;

	COO2ELL(I_precondP,J_precondP,V_precondP,&colELL_precondP, &matrixELL_precondP, 
		&I_COO_LP, &J_COO_LP, &V_COO_LP, numInRowLP, row_idxLP, totalNumPrecond, 
		dimension, &totalNumCOO_LP, maxRowNumPrecondP, &ELL_widthLP);
	printf("ELL_widthLP is %d, and totalNumCOO_LP is %d\n", ELL_widthLP, totalNumCOO_LP);
	cudaMalloc((void **) &col_precondP_d, ELL_widthLP*dimension*sizeof(unsigned int));
	cudaMalloc((void **) &V_precondP_d, ELL_widthLP*dimension*sizeof(double));	
	cudaMemcpy(col_precondP_d,colELL_precondP,dimension*ELL_widthLP*sizeof(unsigned int),cudaMemcpyHostToDevice);
	cudaMemcpy(V_precondP_d,matrixELL_precondP,dimension*ELL_widthLP*sizeof(double),cudaMemcpyHostToDevice);
	if (totalNumCOO_LP>0){
		cudaMalloc((void **) &I_COO_LP_d,totalNumCOO_LP*sizeof(unsigned int));
		cudaMalloc((void **) &J_COO_LP_d,totalNumCOO_LP*sizeof(unsigned int));
		cudaMalloc((void **) &V_COO_LP_d,totalNumCOO_LP*sizeof(double));		
		cudaMemcpy(I_COO_LP_d,I_COO_LP,totalNumCOO_LP*sizeof(unsigned int),cudaMemcpyHostToDevice);
		cudaMemcpy(J_COO_LP_d,J_COO_LP,totalNumCOO_LP*sizeof(unsigned int),cudaMemcpyHostToDevice);
		cudaMemcpy(V_COO_LP_d,V_COO_LP,totalNumCOO_LP*sizeof(double),cudaMemcpyHostToDevice);
	}

	size_t size0=dimension*sizeof(double);
	double *rk_d;//r0 and r1
	double *pk_d;//p0 and p1
	double *bp_d;
	double *zk_d;
	double *zk1_d;
	double *x_d;
	unsigned int* part_boundary_d;
	
	double *rk=(double *)malloc(size0);
	double *zk=(double *)malloc(size0);
	double *zk_1=(double *)malloc(size0);
	
	double *vector_in_d;
	double alphak,betak,dotp0, dotrz0,dotrz1,doth,alphak_1,dotr0;
	
	//double dotz0,dotz0_compare;
	
	double error=10000;
	double threshold;
	cudaMalloc((void **) &rk_d,size0);
	cudaMalloc((void **) &pk_d,size0);
	cudaMalloc((void **) &zk_d,size0);
	cudaMalloc((void **) &zk1_d,size0);
	cudaMalloc((void **) &bp_d,size0);
	cudaMalloc((void **) &x_d,size0);
	cudaMalloc((void **) &part_boundary_d, (partition_size +1)*sizeof(unsigned int));

	cudaMalloc((void **) &vector_in_d,size0);
	//cudaMalloc((void **) &vector_out_d,size0);
	
	cudaMemcpy(vector_in_d, vector_in, size0, cudaMemcpyHostToDevice);	
	if(RODR)
		cudaMemcpy(part_boundary_d, part_boundary, (partition_size +1)*sizeof(unsigned int), cudaMemcpyHostToDevice);	
	
	//
	initialize_bp(dimension,zk_d);
	initialize_bp(dimension,zk1_d);
	initialize_r(dimension, rk_d, vector_in_d);
		

	
	cublasHandle_t handle;
	cublasCreate(&handle);			
	
	if(!BLOCK){
		if(RODR){
			matrix_vectorELL(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
					true, partition_size, part_boundary_d);
		}
		else{
			matrix_vectorELL(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
					false, 0, NULL);
		}
	} else {
		if(RODR){
			matrix_vectorELL_block(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
					true, partition_size, part_boundary_d);
		}
		else{
			matrix_vectorELL_block(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
					false, 0, NULL);
		}
	}

	cudaMemcpy(zk_1, zk1_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
		for(int i = 0; i < 10; ++i)
			printf("zk1[%d] of GPU result is %f\n", i, zk_1[i]);

	if (totalNumCOO_L>0) matrix_vectorCOO(totalNumCOO_L, I_COO_L_d, J_COO_L_d, V_COO_L_d, rk_d, zk1_d, 0, 0);
	if (!BLOCK) {
		if(RODR){
			matrix_vectorELL(dimension, dimension, ELL_widthLP, col_precondP_d, 
					V_precondP_d, zk1_d, zk_d, 0, 0, true, partition_size, part_boundary_d);
		} else {
			matrix_vectorELL(dimension, dimension, ELL_widthLP, col_precondP_d, 
					V_precondP_d, zk1_d, zk_d, 0, 0, false, 0, NULL);
		}
	} else {
		if(RODR){
			matrix_vectorELL_block(dimension, dimension, ELL_widthLP, col_precondP_d, 
					V_precondP_d, zk1_d, zk_d, 0, 0, true, partition_size, part_boundary_d);
		} else {
			matrix_vectorELL_block(dimension, dimension, ELL_widthLP, col_precondP_d, 
					V_precondP_d, zk1_d, zk_d, 0, 0, false, 0, NULL);
		}	
	}
	if (totalNumCOO_LP>0) matrix_vectorCOO(totalNumCOO_LP, I_COO_LP_d, J_COO_LP_d, V_COO_LP_d, zk1_d, zk_d,0,0);
	
	//cublasSdot(handle,dimension,zk_d,1,zk_d,1,&dotz0);
	//for (unsigned int i=0;i<5;i++) printf("dotz0 is %f, dotz0_compare is %f\n",dotz0,dotz0_compare);
	//give an initial value as r0, p0, and set all 
	initialize_all(dimension,pk_d,bp_d,x_d,zk_d,vector_in_d);
	
	
	cublasSdot(handle, dimension,vector_in_d,1,vector_in_d,1,&doth);
	double errorNorm=sqrt(doth);
	unsigned int iter=0;
	struct timeval start1, end1;
	struct timeval start2, end2;

	gettimeofday(&start1, NULL);
	while (iter<MAXIter&&error>1E-8){
		//matrix-vector multiplication, accelerated by our own functions
		if(!BLOCK) {
			if(RODR){
				matrix_vectorELL(dimension, dimension, ELL_width, col_d,V_d,pk_d,bp_d,0,0,
						true,partition_size, part_boundary_d);
			}
			else{
				matrix_vectorELL(dimension, dimension, ELL_width, col_d,V_d,pk_d,bp_d,0,0,
						false, 0, NULL);
			}
		} else {
			if(RODR){
				matrix_vectorELL_block(dimension, dimension, ELL_width, col_d,V_d,pk_d,bp_d,0,0,
						true,partition_size, part_boundary_d);
			}
			else{
				matrix_vectorELL_block(dimension, dimension, ELL_width, col_d,V_d,pk_d,bp_d,0,0,
						false, 0, NULL);
			}
		}
		if (totalNumCOO>0) matrix_vectorCOO(totalNumCOO, I_COO_d, J_COO_d, V_COO_d, pk_d, bp_d, 0, 0);
		//vector-vectror multiplication
		cublasSdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		//r_k*z_k
		cublasSdot(handle,dimension,rk_d,1,zk_d,1,&dotrz0);
		/*cudaMemcpy(rk, rk_d, size0, cudaMemcpyDeviceToHost);
		cudaMemcpy(zk, zk_d, size0, cudaMemcpyDeviceToHost);
		double dotrz0_compare=0;
		for (unsigned int k=0;k<dimension;k++)
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
		if(!BLOCK){
			if(RODR){
				matrix_vectorELL(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
						true, partition_size, part_boundary_d);
			}
			else{
				matrix_vectorELL(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
						false, 0, NULL);
			}
		} else {
			if(RODR){
				matrix_vectorELL_block(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
						true, partition_size, part_boundary_d);
			}
			else{
				matrix_vectorELL_block(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,0,0,
						false, 0, NULL);
			}
		}

		if (totalNumCOO_L>0) matrix_vectorCOO(totalNumCOO_L, I_COO_L_d, J_COO_L_d, V_COO_L_d, rk_d,zk1_d,0 ,0);

		if (!BLOCK){
			if(RODR){	
				matrix_vectorELL(dimension, dimension, ELL_widthLP, col_precondP_d,V_precondP_d,
						zk1_d,zk_d,0,0, true, partition_size, part_boundary_d);
			}
			else{
				matrix_vectorELL(dimension, dimension, ELL_widthLP, col_precondP_d,V_precondP_d,
						zk1_d,zk_d,0,0, false, 0, NULL);
			}
		} else {
			if(RODR){	
				matrix_vectorELL_block(dimension, dimension, ELL_widthLP, col_precondP_d,V_precondP_d,
						zk1_d,zk_d,0,0, true, partition_size, part_boundary_d);
			}
			else{
				matrix_vectorELL_block(dimension, dimension, ELL_widthLP, col_precondP_d,V_precondP_d,
						zk1_d,zk_d,0,0, false, 0, NULL);
			}
		}
		if (totalNumCOO_LP>0) matrix_vectorCOO(totalNumCOO_LP, I_COO_LP_d, J_COO_LP_d, V_COO_LP_d, zk1_d,zk_d,0,0 );
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
		//error_track[iter]=error;
		
		//printf("at iter %d, alpha is %f, beta is %f, dotrz0 is %f, dotrz1 is %f, dotp0 is %f\n", iter, alphak, betak, dotrz0, dotrz1, dotp0);
		iter++;
	}
	printf("max iteration is %d with error %e\n",iter, error);
	cublasDestroy(handle);
	cudaMemcpy(vector_out, x_d, size0, cudaMemcpyDeviceToHost);
	if(totalNumCOO > 0){
		free(I_COO);
		free(J_COO);
		free(V_COO);
	}
	if(totalNumCOO_L > 0){
		free(I_COO_L);
		free(J_COO_L);
		free(V_COO_L);
	}
	if(totalNumCOO_LP > 0){
		free(I_COO_LP);
		free(J_COO_LP);
		free(V_COO_LP);
	}
	free(matrixELL);
	free(colELL);
	free(zk);
	free(zk_1);
	free(rk);
	gettimeofday(&end1, NULL);	
	double timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, time is %f ms, GPU Gflops is %f\n ",iter, timeByMs, (1e-9*(totalNum*4+14*dimension)*1000*iter)/timeByMs);
}

void solverPrecondCPU(const unsigned int procNum, const unsigned int dimension, 
		const unsigned int totalNum, const unsigned int *row_idx, const unsigned int *J, 
		const double *V, const unsigned int totalNumPrecond, const unsigned int *row_idxL, 
		const unsigned int *J_precond, const double *V_precond, const unsigned int totalNumPrecondP,
		const unsigned int *row_idxLP, const unsigned int *J_precondP, const double *V_precondP, 
		const double *vector_in, double *vector_out, const unsigned int MAXIter, unsigned int *realIter){

	size_t size0 =dimension*sizeof(double);
	printf("malloc twice?\n");
	double *zk=(double *)malloc(size0);
	double *zk1=(double *)malloc(size0);
	double *pk=(double *)malloc(size0);
	double *bp=(double *)malloc(size0);
	double *rk=(double *)malloc(size0);
	
	double alphak = 0,betak = 0,dotp0 = 0, dotrz0 = 0,dotrz1 = 0,doth = 0,alphak_1 = 0,dotr0 = 0;

	//double dotz0,dotz0_compare;

	double error=10000;
	double threshold;

	for (unsigned int i=0;i<dimension;i++){
		zk[i]=0;
		zk1[i]=0;
		rk[i]=vector_in[i];
	}

	unsigned int *boundary=(unsigned int *)malloc((procNum+1)*sizeof(unsigned int));
	unsigned int *precond_boundary=(unsigned int *)malloc((procNum+1)*sizeof(unsigned int));
	unsigned int *precondP_boundary=(unsigned int *)malloc((procNum+1)*sizeof(unsigned int));
	boundary[0]=0;
	precond_boundary[0]=0;
	precondP_boundary[0]=0;
	boundary[procNum]=dimension;
	precond_boundary[procNum]=dimension;
	precondP_boundary[procNum]=dimension;

	unsigned int stride, stridePrecond,stridePrecondP;

	stride=ceil((double)totalNum/procNum);
	stridePrecond=ceil((double)totalNumPrecond/procNum);
	stridePrecondP=ceil((double)totalNumPrecondP/procNum);
	unsigned int bias, biasPrecond,biasPrecondP;

	bias=1;
	biasPrecond=1;
	biasPrecondP=1;

	unsigned int a,b;
	
	for (unsigned int i=0;i<dimension;i++){

		if (row_idx[i]>bias*stride){
			boundary[bias]=i;
			bias++;
		}

		if (row_idxL[i]>biasPrecond*stridePrecond){
			precond_boundary[biasPrecond]=i;
			biasPrecond++;
		}		

		if (row_idxLP[i]>biasPrecondP*stridePrecondP){
			precondP_boundary[biasPrecondP]=i;
			biasPrecondP++;
		}

	}


	/*#pragma omp parallel for
	  for (unsigned int i=0;i<totalNumPrecond;i++)
	  zk1[I_precond[i]]+=V_precond[i]*rk[J_precond[i]];*/

	unsigned int rank;
	double tempV;
	/*for(unsigned int i = 0; i < dimension; ++i){
		zk1_compare[i] = 0;
		zk_compare[i] = 0;
	}
	
	for (unsigned int i = 0; i < dimension; ++i){
		tempV=0;
		a=row_idxL[i];
		b=row_idxL[i+1];
		for (unsigned int j=a;j<b;j++){
			tempV+=V_precond[j]*rk[J_precond[j]];
		}
		zk1_compare[i]+=tempV;
		//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
	}*/
	#pragma omp parallel private(rank,tempV, a, b)
	{
		rank=omp_get_thread_num();
		for (unsigned int i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
			tempV=0;
			a=row_idxL[i];
			b=row_idxL[i+1];
			for (unsigned int j=a;j<b;j++){
				tempV+=V_precond[j]*rk[J_precond[j]];
			}
			zk1[i]+=tempV;
			//if(abs(zk1[i] - zk1_compare[i]) > 0.001)
				//printf("at %d parallel computing error with serial %f parallel %f\n", i, zk1_compare[i], zk1[i]);
			//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
		}
	}		

	#pragma omp parallel private(rank,tempV, a, b)
	{
		rank=omp_get_thread_num();
		for (unsigned int i=precondP_boundary[rank];i<precondP_boundary[rank+1];i++){
			tempV=0;
			a=row_idxLP[i];
			b=row_idxLP[i+1];
			for (unsigned int j=a;j<b;j++)
				tempV+=V_precondP[j]*zk1[J_precondP[j]];
			zk[i]+=tempV;
		}

	}				

	//initialize_all(dimension,pk_d,bp_d,x_d,zk_d,vector_in_d);
	#pragma omp parallel for
	for (unsigned int i=0;i<dimension;i++){
		pk[i]=zk[i];
		bp[i]=0;
		vector_out[i]=0;
	}


	//cublasSdot(handle, dimension,vector_in_d,1,vector_in_d,1,&doth);
	#pragma omp parallel for reduction(+:doth)
	for (unsigned int i=0;i<dimension;i++){
		doth+=vector_in[i]*vector_in[i];
	}
	double errorNorm=sqrt(doth);
	//printf("errorNorm is %f\n",errorNorm);
	unsigned int iter=0;
	
	struct timeval start1, end1;
	struct timeval start2, end2;

	gettimeofday(&start1, NULL);
	/*-----------Start the iteration by a while loop-----------*/
	while (iter<MAXIter&&error>1E-8){
		//matrix-vector multiplication, accelerated by our own functions
		/*#pragma omp parallel for
		  for (unsigned int i=0;i<totalNum_1;i++)
		  bp[I_1[i]]+=V_1[i]*pk[J_1[i]];*/
		#pragma omp parallel private(rank,tempV, a, b)
		{
			rank=omp_get_thread_num();
			//gettimeofday(&start2, NULL);
			for (unsigned int i=boundary[rank];i<boundary[rank+1];i++){
				tempV=0;
				a=row_idx[i];
				b=row_idx[i+1];
				for (unsigned int j=a;j<b;j++){
					tempV+=V[j]*pk[J[j]];
				}
				bp[i]+=tempV;
			}
			//gettimeofday(&end2, NULL);
			//printf("at iter %d, rank is %d, boundary span is %d, time is %d us\n", iter, rank, boundary[rank+1]-boundary[rank], ((end2.tv_sec * 1000000 + end2.tv_usec)-(start2.tv_sec * 1000000 + start2.tv_usec)));
		}

		dotp0=0;
		#pragma omp parallel for reduction(+:dotp0)
		for (unsigned int i=0;i<dimension;i++)
			dotp0+=bp[i]*pk[i];

		//r_k*z_k
		dotrz0=0;
		#pragma omp parallel for reduction(+:dotrz0)
		for (unsigned int i=0;i<dimension;i++)
			dotrz0+=rk[i]*zk[i];

		alphak=dotrz0/dotp0;
		alphak_1=-alphak;

		//update x_k to x_k+1, x=x+alphak*p;
		//update r_k to r_k+1, r=r-alphak*bp;
		#pragma omp parallel for
		for (unsigned int i=0;i<dimension;i++){
			vector_out[i]+=alphak*pk[i];
			rk[i]+=alphak_1*bp[i];
			zk1[i]=0;
			zk[i]=0;
		}

		//SpMV: zk=inv(m)*rk
		#pragma omp parallel private(rank,tempV, a, b)
		{
			rank=omp_get_thread_num();
			for (unsigned int i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
				tempV=0;
				a=row_idxL[i];
				b=row_idxL[i+1];
				for (unsigned int j=a;j<b;j++){
					tempV+=V_precond[j]*rk[J_precond[j]];
				}
				zk1[i]+=tempV;
				//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
			}
		}		

		#pragma omp parallel private(rank,tempV, a, b)
		{
			rank=omp_get_thread_num();
			for (unsigned int i=precondP_boundary[rank];i<precondP_boundary[rank+1];i++){
				tempV=0;
				a=row_idxLP[i];
				b=row_idxLP[i+1];
				for (unsigned int j=a;j<b;j++)
					tempV+=V_precondP[j]*zk1[J_precondP[j]];
				zk[i]+=tempV;
			}

		}			

		//r_k+1 * z_k+1
		dotrz1=0;
		#pragma omp parallel for reduction(+:dotrz1)
		for (unsigned int i=0;i<dimension;i++)
			dotrz1+=rk[i]*zk[i];
		betak=dotrz1/dotrz0;
		//printf("dotp0 is %f, dotrz0 is %f dotrz1 is %f, betak is %f and alphak is %f at iter %d\n", dotp0, dotrz0,dotrz1,betak, alphak, iter);
		//p=r+gamak*p;
		#pragma omp parallel for
		for (unsigned int i=0;i<dimension;i++){
			pk[i]=zk[i]+betak*pk[i];
			bp[i]=0;
		}
		dotr0=0;
		#pragma omp parallel for reduction (+:dotr0)
		for (unsigned int i=0;i<dimension;i++)
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
	double timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, Xeon_Phi Gflops is %f\n ",iter, (1e-9*(totalNum*4+14*dimension)*1000*iter)/timeByMs);



}


