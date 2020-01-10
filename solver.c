#include "kernel.h"
#include "solver.h"
#include "convert.h"
#include "test.h"

static inline uint32_t calPadding(uint32_t dimension,
		const uint32_t* numInRow, const uint32_t* part_boundary, 	
		uint32_t* ELL_block_cols_vec, const uint32_t part_size){

	uint32_t padded = 0;
	for(uint32_t i = 0; i < part_size; ++i){
		for(uint32_t j = part_boundary[i]; j < part_boundary[i+1]; ++j){
			if(j < part_boundary[i] + block_per_part*ELL_threadSize){
				uint32_t blockIdx = i*block_per_part + (j-part_boundary[i])/ELL_threadSize;
				if(numInRow[j] < ELL_block_cols_vec[blockIdx])	
					padded += (ELL_block_cols_vec[blockIdx] - numInRow[j]);
			}
		}
	}
	return padded;
}


static void ELL_cuda_malloc_trans_data(uint32_t** col_d, double** V_d, 
				uint32_t* colELL, double* matrixELL,
				const uint32_t dimension, const uint32_t ELL_width){

	cudaMalloc((void **) col_d, ELL_width*dimension*sizeof(uint32_t));
    cudaMalloc((void **) V_d, ELL_width*dimension*sizeof(double));
    cudaMemcpy(*col_d,colELL,dimension*ELL_width*sizeof(uint32_t),cudaMemcpyHostToDevice);
    cudaMemcpy(*V_d,matrixELL,dimension*ELL_width*sizeof(double),cudaMemcpyHostToDevice);
}

static void ELL_cuda_malloc_trans_data_block(uint32_t** col_d, double** V_d, 
				uint32_t** ELL_block_cols_vec_d, uint32_t** ELL_block_bias_vec_d,
                uint32_t* colELL, double* matrixELL,
				uint32_t* ELL_block_cols_vec, uint32_t* ELL_block_bias_vec,
				uint32_t dimension, uint32_t blocks){
                
	uint32_t ELL_size = ELL_block_bias_vec[blocks];
    cudaMalloc((void **) col_d, ELL_size*sizeof(uint32_t));
    cudaMalloc((void **) V_d, ELL_size*sizeof(double));
    cudaMalloc((void **) ELL_block_cols_vec_d, blocks*sizeof(uint32_t));
    cudaMalloc((void **) ELL_block_bias_vec_d, (blocks + 1)*sizeof(uint32_t));

    cudaMemcpy(*col_d,colELL,ELL_size*sizeof(uint32_t),cudaMemcpyHostToDevice);
    cudaMemcpy(*V_d,matrixELL,ELL_size*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(*ELL_block_cols_vec_d, ELL_block_cols_vec, blocks*sizeof(uint32_t),cudaMemcpyHostToDevice);
    cudaMemcpy(*ELL_block_bias_vec_d, ELL_block_bias_vec, (blocks + 1)*sizeof(uint32_t),cudaMemcpyHostToDevice);
}



static void COO_cuda_malloc_trans_data(uint32_t** I_COO_d, uint32_t** J_COO_d, double** V_COO_d,
				uint32_t*  I_COO, uint32_t* J_COO, double* V_COO,
				const uint32_t dimension, const uint32_t totalNumCOO)
{
	cudaMalloc((void **) I_COO_d,totalNumCOO*sizeof(uint32_t));
	cudaMalloc((void **) J_COO_d,totalNumCOO*sizeof(uint32_t));
	cudaMalloc((void **) V_COO_d,totalNumCOO*sizeof(double));	
	cudaMemcpy(*I_COO_d,I_COO,totalNumCOO*sizeof(uint32_t),cudaMemcpyHostToDevice);
	cudaMemcpy(*J_COO_d,J_COO,totalNumCOO*sizeof(uint32_t),cudaMemcpyHostToDevice);
	cudaMemcpy(*V_COO_d,V_COO,totalNumCOO*sizeof(double),cudaMemcpyHostToDevice);
}

void solver(const uint32_t dimension, const uint32_t totalNum, const uint32_t *I, 
			const uint32_t *J, const double *V, const double *vector_in, 
			double *vector_out, double *error_track, uint32_t MAXIter)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for
	double dotp0,dotr0,dotr1,doth;
	size_t size1=dimension*sizeof(double);
	double *bp=(double *) malloc(size1);
	double *pk=(double *) malloc(size1);
	double *rk=(double *) malloc(size1);
	//double *x=(double *) malloc(size1);
	uint32_t i;
	double threshold=0.0000001;
	uint32_t iter=0;
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
			//for (uint32_t k=0;k<5;k++) printf("pk %d at iter %d is %f\n",k, iter, pk[k]);
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

void solverGPU_unprecondHYB(matrixCOO_S* localMatrix, 
		const double *vector_in, double *vector_out,  
		const uint32_t MAXIter, uint32_t *realIter,  const cb_s cb,
		const uint32_t part_size, const uint32_t* part_boundary)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for
	double dotp0,dotr0,dotr1,doth;
	uint32_t dimension, totalNum, maxRowNum, totalNumCOO; 
    uint32_t *row_idx, *numInRow, *I, *J; 
    double *V;
    dimension = localMatrix->dimension; 
    totalNum = localMatrix->totalNum; 
    maxRowNum = localMatrix->maxRowNum; 

	uint32_t ELL_blocks;
    uint32_t *ELL_block_cols_vec, *ELL_block_cols_vec_d;    
    uint32_t *ELL_block_bias_vec, *ELL_block_bias_vec_d;	
	ELL_block_cols_vec = (uint32_t*)malloc(ELL_blocks*sizeof(uint32_t));  
	ELL_block_bias_vec = (uint32_t*)malloc((ELL_blocks +1)*sizeof(uint32_t));

    row_idx = localMatrix->row_idx; 
    numInRow = localMatrix->numInRow; 
    I = localMatrix->I; 
    J = localMatrix->J;
    V = localMatrix->V;
	uint32_t *colELL, *I_COO, *J_COO, ELL_width;
	double* matrixELL, *V_COO;
	uint32_t *col_d;
	double *V_d;
	volatile bool RODR, BLOCK, CACHE;
	uint32_t *I_COO_d, *J_COO_d;
	double *V_COO_d;
	BLOCK = cb.BLOCK;
	CACHE = cb.CACHE;
	RODR = cb.RODR;
	uint32_t *part_boundary_d;
	
	if(BLOCK){
        //RODR will change the number of blocks since it provides cache option
        if(RODR)
			ELL_blocks = part_size*block_per_part; 
		else
			ELL_blocks = ceil((double) dimension/ELL_threadSize);

        ELL_block_cols_vec = (uint32_t*)malloc(ELL_blocks*sizeof(uint32_t));  
        ELL_block_bias_vec = (uint32_t*)malloc((ELL_blocks +1)*sizeof(uint32_t));  
		COO2ELL_block(&totalNumCOO, 
				ELL_block_cols_vec, 
				ELL_block_bias_vec,
				&colELL, 
				&matrixELL, 
				&I_COO, 
				&J_COO, 
				&V_COO,
				I, 
				J, 
				V, 
				row_idx, 
				numInRow, 
				maxRowNum, 
				totalNum, 
				dimension,
				part_size, 
				ELL_blocks,
				part_boundary, 
				RODR);
		ELL_cuda_malloc_trans_data_block(&col_d, 
				&V_d, 
				&ELL_block_cols_vec_d, 
				&ELL_block_bias_vec_d,
				colELL, 
				matrixELL,
				ELL_block_cols_vec, 
				ELL_block_bias_vec,
				dimension, 
				ELL_blocks);

		cudaMalloc((void **) &part_boundary_d, (part_size +1)*sizeof(uint32_t));
		cudaMemcpy(part_boundary_d, part_boundary, 
				(part_size +1)*sizeof(uint32_t), cudaMemcpyHostToDevice);
		uint32_t padding = calPadding(dimension, numInRow, part_boundary, 	
			ELL_block_cols_vec, part_size);
		printf("padding is %d, waste rate on ELL is %f\n", padding, ((float)padding)/totalNum);

	} else {
		COO2ELL(I,J,V,&colELL,&matrixELL,&I_COO, &J_COO, &V_COO, 
			numInRow, row_idx, totalNum, dimension, &totalNumCOO, maxRowNum, &ELL_width);

        ELL_cuda_malloc_trans_data(&col_d, &V_d, colELL, matrixELL, dimension, ELL_width);
	}
	
	if (totalNumCOO>0){
        COO_cuda_malloc_trans_data(&I_COO_d, &J_COO_d, &V_COO_d, 
				I_COO, J_COO, V_COO, 
				dimension, totalNumCOO);
	}
	cublasHandle_t handle;
	cublasCreate(&handle);

	matrixHYB_S_d localMatrixHYB_d;
	
	localMatrixHYB_d.dimension = dimension;
	localMatrixHYB_d.ELL_width = ELL_width;
	localMatrixHYB_d.totalNumCOO = totalNumCOO;
	localMatrixHYB_d.col_d = col_d;
	localMatrixHYB_d.I_COO_d = I_COO_d;
	localMatrixHYB_d.J_COO_d = J_COO_d;
	localMatrixHYB_d.V_d = V_d;
	localMatrixHYB_d.V_COO_d = V_COO_d;
	localMatrixHYB_d.ELL_block_bias_vec_d = ELL_block_bias_vec_d;
	localMatrixHYB_d.ELL_block_cols_vec_d = ELL_block_cols_vec_d;

	double *bp_d, *pk_d, *rk_d, *vector_out_d;
	size_t size0=dimension*sizeof(uint32_t);
	size_t size1=dimension*sizeof(double);
	cudaMalloc((void **) &bp_d,size1);
	cudaMalloc((void **) &pk_d,size1);
	cudaMalloc((void **) &rk_d,size1);
	cudaMalloc((void **) &vector_out_d,size1);
	//double *x=(double *) malloc(size1);
	double threshold=0.0000001;
	uint32_t iter=0;
	double const1 = 1.0;
	double error, alphak, _alphak, gamak;
	error=1000;
	//initialize
	doth=0;
    cudaMemcpy(pk_d, vector_in, dimension*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(rk_d, vector_in, dimension*sizeof(double), cudaMemcpyHostToDevice);
	for (uint32_t i=0;i<dimension;i++) {
		doth=doth+vector_in[i]*vector_in[i];
	}
	struct timeval start1, end1;
	gettimeofday(&start1, NULL);
	double *bp=(double *) malloc(size1);
	double *pk=(double *) malloc(size1);
	double *rk=(double *) malloc(size1);

	double *bp_dt = (double *) malloc(size1);
	double *pk_dt = (double *) malloc(size1);
	double *rk_dt = (double *) malloc(size1);
	//double *x=(double *) malloc(size1);
	double error_t,alphak_t,gamak_t;
	error=1000;
	//initialize
	for (uint32_t i=0;i<dimension;i++)
	{
		pk[i]=vector_in[i];
		rk[i]=vector_in[i];
		vector_out[i]=0;
		bp[i]=0;
	}
	
	while (error>threshold&&iter<MAXIter){
		dotp0=0;
		dotr0=0;
		dotr1=0;
		uint32_t errorIdx = 0;
		double compareError;
		cudaMemset(bp_d, 0, size1);
		matrix_vectorHYB(&localMatrixHYB_d, pk_d, bp_d, cb, 0,
			   part_size, part_boundary_d);
		
		//for (uint32_t k=0;k<5;k++) printf("pk %d at iter %d is %f\n",k, iter, pk[k]);
		//printf("CPU start\n");
		
		cublasDdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr0);
			
		alphak=dotr0/dotp0;
		_alphak = -alphak;
		
		cublasDaxpy(handle,dimension,&alphak,pk_d,1,vector_out_d,1);
		cublasDaxpy(handle,dimension,&_alphak,bp_d,1,rk_d,1);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr1);
		
		gamak=dotr1/dotr0;

		cublasDscal(handle,dimension,&gamak,pk_d,1);
		cublasDaxpy(handle,dimension,&const1, rk_d, 1, pk_d, 1);
		
		//printf("at iter %d, alphak is %f, gamak is %f\n",iter, alphak,gamak);
		error=sqrt(dotr1)/sqrt(doth);
		//error_track[iter]=error;
		//printf("error at %d is %f\n",iter, error);
		iter++;
	}
	cudaMemcpy(vector_out, vector_out_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
	gettimeofday(&end1, NULL);	
	double timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, time is %f ms, GPU Gflops is %f\n ",iter, timeByMs, 
			(1e-9*(totalNum*2+13*dimension)*1000*iter)/timeByMs);
}

void solverGPU_unprecondCUSPARSE(matrixCOO_S* localMatrix, 
		const double *vector_in, double *vector_out,  
		const uint32_t MAXIter, uint32_t *realIter,  const cb_s cb,
		const uint32_t part_size, const uint32_t* part_boundary)
{
	//exampine the performance using cusparse library functions with
	//CSR format
	double dotp0,dotr0,dotr1,doth;

}




static uint32_t get_blocks_rodr(uint32_t* part_boundary, const uint32_t part_size){
	uint32_t block_num = 0;
	for(uint32_t i = 0; i < part_size; ++i){
		ceil(((float) shared_per_block)/(ELL_threadSize*element_size));
		block_num += ceil((float (part_boundary[i + 1] - part_boundary[i]))/ELL_threadSize);		
	}
	return block_num;
}

void solverGPU_HYB(matrixCOO_S* localMatrix, matrixCOO_S* localMatrix_precond, 
                matrixCOO_S* localMatrix_precondP,
        		const double *vector_in, double *vector_out,  
                const uint32_t MAXIter, uint32_t *realIter,  const cb_s cb,
                const uint32_t part_size, const uint32_t* part_boundary)
{
    uint32_t dimension, totalNum, maxRowNum, totalNumPrecond; 
    uint32_t maxRowNumPrecond, totalNumPrecondP, maxRowNumPrecondP;
    uint32_t *row_idx, *numInRow, *I, *J; 
    double *V;
    uint32_t *numInRowL, *row_idxL; 
    const uint32_t *I_precond, *J_precond; 
    const double* V_precond; 
    const uint32_t *numInRowLP; 
    const uint32_t *row_idxLP;  
    const uint32_t *I_precondP, *J_precondP; 
    const double*V_precondP; 
	bool FACT = cb.FACT;
    dimension = localMatrix->dimension; 
    totalNum = localMatrix->totalNum; 
    maxRowNum = localMatrix->maxRowNum; 

    totalNumPrecond = localMatrix_precond->totalNum;
    maxRowNumPrecond = localMatrix_precond->maxRowNum; 
	if(FACT){	
		totalNumPrecondP = localMatrix_precondP->totalNum; 
		maxRowNumPrecondP = localMatrix_precondP->maxRowNum;
	}

    row_idx = localMatrix->row_idx; 
    numInRow = localMatrix->numInRow; 
    I = localMatrix->I; 
    J = localMatrix->J;
    V = localMatrix->V;

    row_idxL = localMatrix_precond->row_idx; 
    numInRowL = localMatrix_precond->numInRow; 
    I_precond = localMatrix_precond->I; 
    J_precond = localMatrix_precond->J;
    V_precond = localMatrix_precond->V;
	if(FACT){
		row_idxLP = localMatrix_precondP->row_idx; 
		numInRowLP = localMatrix_precondP->numInRow; 
		I_precondP = localMatrix_precondP->I; 
		J_precondP = localMatrix_precondP->J;
		V_precondP = localMatrix_precondP->V;
	}


	uint32_t *colELL, *I_COO, *J_COO, ELL_width, ELL_widthL, ELL_widthLP;
	double* matrixELL, *V_COO;
	uint32_t *col_d;
	double *V_d;
	volatile bool RODR, BLOCK, CACHE;
	uint32_t *I_COO_d, *J_COO_d;
	double *V_COO_d;
	BLOCK = cb.BLOCK;
	CACHE = cb.CACHE;
	RODR = cb.RODR;
    uint32_t ELL_blocks;
    uint32_t *ELL_block_cols_vec, *ELL_block_cols_vec_d;    
    uint32_t *ELL_block_bias_vec, *ELL_block_bias_vec_d;	
	uint32_t *ELL_block_cols_vec_L, *ELL_block_cols_vec_L_d;    
    uint32_t *ELL_block_bias_vec_L, *ELL_block_bias_vec_L_d;
	uint32_t *ELL_block_cols_vec_LP, *ELL_block_cols_vec_LP_d;    
	uint32_t *ELL_block_bias_vec_LP, *ELL_block_bias_vec_LP_d;
	uint32_t totalNumCOO, totalNumCOO_L, totalNumCOO_LP;
    if(BLOCK){
        //RODR will change the number of blocks since it provides cache option
        if(RODR)
			ELL_blocks = part_size*block_per_part; 
		else
			ELL_blocks = ceil((double) dimension/ELL_threadSize);

        ELL_block_cols_vec = (uint32_t*)malloc(ELL_blocks*sizeof(uint32_t));  
        ELL_block_bias_vec = (uint32_t*)malloc((ELL_blocks +1)*sizeof(uint32_t));  
		COO2ELL_block(&totalNumCOO, 
				ELL_block_cols_vec, 
				ELL_block_bias_vec,
				&colELL, 
				&matrixELL, 
				&I_COO, 
				&J_COO, 
				&V_COO,
				I, 
				J, 
				V, 
				row_idx, 
				numInRow, 
				maxRowNum, 
				totalNum, 
				dimension,
				part_size, 
				ELL_blocks,
				part_boundary, 
				RODR);
		ELL_cuda_malloc_trans_data_block(&col_d, 
				&V_d, 
				&ELL_block_cols_vec_d, 
				&ELL_block_bias_vec_d,
				colELL, 
				matrixELL,
				ELL_block_cols_vec, 
				ELL_block_bias_vec,
				dimension, 
				ELL_blocks);
	} else {
		COO2ELL(I,J,V,&colELL,&matrixELL,&I_COO, &J_COO, &V_COO, 
			numInRow, row_idx, totalNum, dimension, &totalNumCOO, maxRowNum, &ELL_width);

        ELL_cuda_malloc_trans_data(&col_d, &V_d, colELL, matrixELL, dimension, ELL_width);
	}
	
	if (totalNumCOO>0){
        COO_cuda_malloc_trans_data(&I_COO_d, &J_COO_d, &V_COO_d, 
				I_COO, J_COO, V_COO, 
				dimension, totalNumCOO);
	}	

	/*matrix L'*/
	uint32_t *colELL_precond, *I_COO_L, *J_COO_L;
	double* matrixELL_precond, *V_COO_L;
	uint32_t *col_precond_d;
	double *V_precond_d;
	uint32_t *I_COO_L_d, *J_COO_L_d;
	double *V_COO_L_d;

    if(BLOCK){
		ELL_block_cols_vec_L = (uint32_t*)malloc(ELL_blocks*sizeof(uint32_t));  
        ELL_block_bias_vec_L = (uint32_t*)malloc((ELL_blocks + 1)*sizeof(uint32_t));
		COO2ELL_block(&totalNumCOO_L, 
				ELL_block_cols_vec_L, 
				ELL_block_bias_vec_L,
				&colELL_precond, 
				&matrixELL_precond,
			   	&I_COO_L, &J_COO_L, &V_COO_L,
				I_precond, J_precond, V_precond, 
				row_idxL, numInRowL, maxRowNumPrecond,
				totalNumPrecond, dimension, part_size,
				ELL_blocks, part_boundary, RODR);
		ELL_cuda_malloc_trans_data_block(&col_precond_d, &V_precond_d, 
				&ELL_block_cols_vec_L_d, &ELL_block_bias_vec_L_d,
				colELL_precond, matrixELL_precond,
				ELL_block_cols_vec_L, ELL_block_bias_vec_L,
				dimension, ELL_blocks);
	} else {
		COO2ELL(I_precond,J_precond,V_precond,&colELL_precond, &matrixELL_precond, 
				&I_COO_L, &J_COO_L, &V_COO_L, numInRowL, row_idxL, totalNumPrecond, 
				dimension, &totalNumCOO_L, maxRowNumPrecond, &ELL_widthL);
        ELL_cuda_malloc_trans_data(&col_precond_d, &V_precond_d, colELL_precond, 
				matrixELL_precond, dimension, ELL_widthL);
	}
	//printf("ELL_widthL is %d, and totalNumCOO_L is %d\n", ELL_widthL, totalNumCOO_L);
	if (totalNumCOO_L>0){
		COO_cuda_malloc_trans_data(&I_COO_L_d, &J_COO_L_d, &V_COO_L_d, 
				I_COO_L, J_COO_L, V_COO_L, 
				dimension, totalNumCOO_L);
	}	

	/*matrix L*/
	uint32_t *colELL_precondP, *I_COO_LP, *J_COO_LP;
	double* matrixELL_precondP, *V_COO_LP;
	uint32_t *col_precondP_d;
	double *V_precondP_d;
	uint32_t *I_COO_LP_d, *J_COO_LP_d;
	double *V_COO_LP_d;
	if(FACT){
		if(BLOCK){
			ELL_block_cols_vec_LP = (uint32_t*)malloc(ELL_blocks*sizeof(uint32_t));  
			ELL_block_bias_vec_LP = (uint32_t*)malloc((ELL_blocks + 1)*sizeof(uint32_t));
			COO2ELL_block(&totalNumCOO_LP, ELL_block_cols_vec_LP, ELL_block_bias_vec_LP,
					&colELL_precondP, &matrixELL_precondP, &I_COO_LP, &J_COO_LP, &V_COO_LP,
					I_precondP, J_precondP, V_precondP, 
					row_idxLP, numInRowLP, maxRowNumPrecondP, 
					totalNumPrecondP, dimension, part_size,
					ELL_blocks, part_boundary, RODR);
			printf("the last cols_vec_LP is %d\n", ELL_block_cols_vec_LP[ELL_blocks - 1]);
			ELL_cuda_malloc_trans_data_block(&col_precondP_d, &V_precondP_d, 
					&ELL_block_cols_vec_LP_d, &ELL_block_bias_vec_LP_d,
					colELL_precondP, matrixELL_precondP,
					ELL_block_cols_vec_LP, ELL_block_bias_vec_LP,
					dimension, ELL_blocks);
		} else {
			COO2ELL(I_precondP,J_precondP,V_precondP,&colELL_precondP, &matrixELL_precondP, 
					&I_COO_LP, &J_COO_LP, &V_COO_LP, numInRowLP, row_idxLP, totalNumPrecondP, 
					dimension, &totalNumCOO_LP, maxRowNumPrecondP, &ELL_widthLP);
			ELL_cuda_malloc_trans_data(&col_precondP_d, &V_precondP_d, colELL_precondP, 
					matrixELL_precondP, dimension, ELL_widthLP);
		}
		//printf("ELL_widthL is %d, and totalNumCOO_L is %d\n", ELL_widthL, totalNumCOO_L);
		if (totalNumCOO_LP>0){
			COO_cuda_malloc_trans_data(&I_COO_LP_d, &J_COO_LP_d, &V_COO_LP_d, 
					I_COO_LP, J_COO_LP, V_COO_LP, 
					dimension, totalNumCOO_LP);
		}
	}	
	printf("totalNumCOO is  %d totalNumCOO_L is %d, and totalNumCOO_LP is %d\n", totalNumCOO, totalNumCOO_L, totalNumCOO_LP);
	//printf("fill rate original is %f, fill rate L is %f, fill rate lp is %f\n", calFilledRate(totalNum, ELL_blocks,), calFilledRate, calFilledRate);
	
	size_t size0=dimension*sizeof(double);
	double *rk_d;//r0 and r1
	double *pk_d;//p0 and p1
	double *bp_d;
	double *zk_d;
	double *zk1_d;
	double *x_d;
	uint32_t* part_boundary_d;
	
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
	if(BLOCK && RODR){
	   	cudaMalloc((void **) &part_boundary_d, (part_size +1)*sizeof(uint32_t));
		cudaMemcpy(part_boundary_d, part_boundary, 
				(part_size +1)*sizeof(uint32_t), cudaMemcpyHostToDevice);
	}

	cudaMalloc((void **) &vector_in_d,size0);
	//cudaMalloc((void **) &vector_out_d,size0);
	cudaMemcpy(vector_in_d, vector_in, size0, cudaMemcpyHostToDevice);	
	
	//
	initialize_bp(dimension,zk_d);
	initialize_bp(dimension,zk1_d);
	initialize_r(dimension, rk_d, vector_in_d);

	double* pkT = (double*) malloc(dimension*sizeof(double));	
	double* bpT = (double*) calloc(dimension, sizeof(double));	
	double* bp_g = (double*) calloc(dimension, sizeof(double));

	cublasHandle_t handle;
	cublasCreate(&handle);			
	double* zk1T = (double*) calloc(dimension, sizeof(double));	
	double* rkT = (double*) calloc(dimension, sizeof(double));	
	double* zk1_g = (double*) calloc(dimension, sizeof(double));
	if(!BLOCK){
		matrix_vectorELL(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,
				rk_d,zk1_d, false, 0, NULL);
					
	} else {
		if(RODR){
			matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_L_d, 
					ELL_block_bias_vec_L_d,
                    col_precond_d,V_precond_d,rk_d,zk1_d,
					CACHE, part_size, part_boundary_d);
		}
		else{
			matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_L_d, 
					ELL_block_bias_vec_L_d,
					col_precond_d,V_precond_d,rk_d,zk1_d,
					false, 0, NULL);
		}
	}

	if (totalNumCOO_L>0) matrix_vectorCOO(totalNumCOO_L, I_COO_L_d, J_COO_L_d, V_COO_L_d, 
									rk_d, zk1_d);
	cudaMemcpy(zk1_g, zk1_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
	uint32_t errorIdx = 0;
	float compareError;
	matrix_vectorTest(localMatrix_precond, vector_in, zk1T, 0);
	for(int i = 0; i < dimension; ++i){
		rkT[i] = vector_in[i];
		compareError = 	(zk1T[i] - zk1_g[i])/zk1_g[i];
		if(errorIdx == 0  && (compareError > 0.000001 || compareError < -0.00001)){ 	
			printf("zk1[%d] of GPU result is %f, test value is %f\n", i, zk1_g[i], zk1T[i]);
			if(errorIdx == 0) errorIdx = i;	
		}
	}
	if(errorIdx > 0){
		matrix_vectorTest(localMatrix_precond, vector_in, zk1T, errorIdx);
		matrix_vectorELL_block(dimension, errorIdx, ELL_block_cols_vec_L_d, 
				ELL_block_bias_vec_L_d,
				col_precond_d,V_precond_d,rk_d,zk1_d,
				CACHE, part_size, part_boundary_d);
	}
	if(FACT){	
		if (!BLOCK) {
			matrix_vectorELL(dimension, dimension, ELL_widthLP, col_precondP_d, 
					V_precondP_d, zk1_d, zk_d, false, 0, NULL);
		} else {
			if(RODR){
				matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_LP_d, ELL_block_bias_vec_LP_d,  
						col_precondP_d,V_precondP_d, 
						zk1_d, zk_d, CACHE, part_size, part_boundary_d); 

			} else {
				matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_LP_d, ELL_block_bias_vec_LP_d,  
						col_precondP_d, V_precondP_d, 
						zk1_d, zk_d, false, 0, NULL);

			}	
		}
		if (totalNumCOO_LP>0) matrix_vectorCOO(totalNumCOO_LP, I_COO_LP_d, J_COO_LP_d, V_COO_LP_d, zk1_d, zk_d);
	}	
	//double* zkT = (double*) calloc(dimension, sizeof(double));	
	//double* zk_g = (double*) calloc(dimension, sizeof(double));	
	//cudaMemcpy(zk_g, zk_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
	//matrix_vectorTest(localMatrix_precondP, zk1T, zkT);
	//for(int i = 0; i < dimension; ++i){
	//	pkT[i] = zkT[i];
	//	if(abs((zkT[i] - zk_g[i])/zk_g[i]) > 0.000001) 	
	//		printf("zk[%d] of GPU result is %f, test value is %f\n", i, zk_g[i], zkT[i]);
	//}
	//cublasDdot(handle,dimension,zk_d,1,zk_d,1,&dotz0);
	//for (uint32_t i=0;i<5;i++) printf("dotz0 is %f, dotz0_compare is %f\n",dotz0,dotz0_compare);
	//give an initial value as r0, p0, and set all 
	initialize_all(dimension,pk_d,bp_d,x_d,zk_d,vector_in_d);
	
	
	cublasDdot(handle, dimension,vector_in_d,1,vector_in_d,1,&doth);
	double errorNorm=sqrt(doth);
	uint32_t iter=0;
	struct timeval start1, end1;
	struct timeval start2, end2;

	gettimeofday(&start1, NULL);
	

	while (iter<MAXIter&&error>1E-8){
		//matrix-vector multiplication, accelerated by our own functions
		if(!BLOCK) {
			matrix_vectorELL(dimension, dimension, ELL_width, col_d,V_d,pk_d,bp_d,
								false, 0, NULL);
		} else {
			if(RODR){
				matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_d, ELL_block_bias_vec_d,  
						col_d,V_d,pk_d,bp_d,
						CACHE, part_size, part_boundary_d);
			}
			else{
				matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_d, ELL_block_bias_vec_d,  
						col_d,V_d,pk_d,bp_d,
						false, 0, NULL);
			}
		}
		if (totalNumCOO>0) matrix_vectorCOO(totalNumCOO, I_COO_d, J_COO_d, V_COO_d, pk_d, bp_d);
			
		//cudaMemcpy(bp_g, bp_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
		//double dotp0T = 0;
		//double dotp0TT = 0;
		//double dotrz0T = 0;
		//matrix_vectorTest(localMatrix, pkT, bpT);
		//for(int i = 0; i < dimension; ++i){
		//	dotp0T += pkT[i]*bpT[i];
		//	dotp0TT += pkT[i]*bp_g[i];
		//	dotrz0T += rkT[i]*zkT[i];
		//	if(abs(dotp0TT - dotp0T) > 0.000001) 	
		//		printf("dotp0TT[%d] of GPU result is %f, test value is %f\n", i, dotp0TT, dotp0T);
		//}
		//vector-vectror multiplication
		cublasDdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		//r_k*z_k
		cublasDdot(handle,dimension,rk_d,1,zk_d,1,&dotrz0);
		/*cudaMemcpy(rk, rk_d, size0, cudaMemcpyDeviceToHost);
		cudaMemcpy(zk, zk_d, size0, cudaMemcpyDeviceToHost);
		double dotrz0_compare=0;
		for (uint32_t k=0;k<dimension;k++)
		{
			dotrz0_compare=dotrz0_compare+rk[k]*zk[k];
		}		
		printf("dotrz0 is %f, dotrz0_compare is %f\n",dotrz0,dotrz0_compare);*/
		alphak=dotrz0/dotp0;
		//double alphakT = dotrz0T/dotp0T;
		alphak_1=-alphak;
		
		//update x_k to x_k+1, x=x+alphak*p;
		cublasDaxpy(handle,dimension,&alphak,pk_d,1,x_d,1);
		//update r_k to r_k+1, r=r-alphak*p;
		cublasDaxpy(handle,dimension,&alphak_1,bp_d,1,rk_d,1);
		initialize_bp(dimension,zk_d);
		initialize_bp(dimension,zk1_d);
		if(!BLOCK){
			matrix_vectorELL(dimension, dimension, ELL_widthL, col_precond_d,V_precond_d,rk_d,zk1_d,
					false, 0, NULL);
		} else {
			if(RODR){
				matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_L_d, ELL_block_bias_vec_L_d, 
						col_precond_d, V_precond_d, rk_d, zk1_d, 
						CACHE, part_size, part_boundary_d);
			}
			else{
				matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_L_d, ELL_block_bias_vec_L_d, 
						col_precond_d,V_precond_d,rk_d,zk1_d,
						false, 0, NULL);
			}
		}

		if (totalNumCOO_L>0) matrix_vectorCOO(totalNumCOO_L, I_COO_L_d, J_COO_L_d, V_COO_L_d, rk_d,zk1_d);
		if(FACT){
			if (!BLOCK){
				matrix_vectorELL(dimension, dimension, ELL_widthLP, col_precondP_d,V_precondP_d,
						zk1_d,zk_d, false, 0, NULL);
			} else {
				if(RODR){	
					matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_LP_d, ELL_block_bias_vec_LP_d,  
							col_precondP_d,V_precondP_d,
							zk1_d,zk_d, CACHE, part_size, part_boundary_d);
				}
				else{
					matrix_vectorELL_block(dimension, 0, ELL_block_cols_vec_LP_d, ELL_block_bias_vec_LP_d,  
							col_precondP_d,V_precondP_d,
							zk1_d,zk_d, false, 0, NULL);
				}
			}
			if (totalNumCOO_LP>0) matrix_vectorCOO(totalNumCOO_LP, I_COO_LP_d, J_COO_LP_d, V_COO_LP_d, zk1_d,zk_d );
		}
		//r_k+1 * z_k+1
		cublasDdot(handle,dimension,rk_d,1,zk_d,1,&dotrz1);
		betak=dotrz1/dotrz0;
		//printf("dotp0 is %f, dotrz0 is %f dotrz1 is %f, betak is %f and alphak is %f at iter %d\n", dotp0, dotrz0,dotrz1,betak, alphak, iter);
		//p=r+gamak*p;
		myxpy(dimension,betak,zk_d,pk_d);//if using cublas, we need two functin to complete it
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr0);
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
    if(BLOCK){
        free(ELL_block_cols_vec);
        free(ELL_block_bias_vec);
        free(ELL_block_cols_vec_L);
        free(ELL_block_bias_vec_L);
		if(FACT){
			free(ELL_block_cols_vec_LP);
			free(ELL_block_bias_vec_LP);
		}
    }

	gettimeofday(&end1, NULL);	
	double timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, time is %f ms, GPU Gflops is %f\n ",iter, timeByMs, (1e-9*(totalNum*4+14*dimension)*1000*iter)/timeByMs);
}

void solverPrecondCPU(const uint32_t procNum, const uint32_t dimension, 
		const uint32_t totalNum, const uint32_t *row_idx, const uint32_t *J, 
		const double *V, const uint32_t totalNumPrecond, const uint32_t *row_idxL, 
		const uint32_t *J_precond, const double *V_precond, const uint32_t totalNumPrecondP,
		const uint32_t *row_idxLP, const uint32_t *J_precondP, const double* V_precondP, 
		const double *vector_in, double *vector_out, const uint32_t MAXIter, uint32_t *realIter){

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

	for (uint32_t i=0;i<dimension;i++){
		zk[i]=0;
		zk1[i]=0;
		rk[i]=vector_in[i];
	}

	uint32_t *boundary=(uint32_t *)malloc((procNum+1)*sizeof(uint32_t));
	uint32_t *precond_boundary=(uint32_t *)malloc((procNum+1)*sizeof(uint32_t));
	uint32_t *precondP_boundary=(uint32_t *)malloc((procNum+1)*sizeof(uint32_t));
	boundary[0]=0;
	precond_boundary[0]=0;
	precondP_boundary[0]=0;
	boundary[procNum]=dimension;
	precond_boundary[procNum]=dimension;
	precondP_boundary[procNum]=dimension;

	uint32_t stride, stridePrecond,stridePrecondP;

	stride=ceil((double)totalNum/procNum);
	stridePrecond=ceil((double)totalNumPrecond/procNum);
	stridePrecondP=ceil((double)totalNumPrecondP/procNum);
	uint32_t bias, biasPrecond,biasPrecondP;

	bias=1;
	biasPrecond=1;
	biasPrecondP=1;

	uint32_t a,b;
	
	for (uint32_t i=0;i<dimension;i++){

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
	  for (uint32_t i=0;i<totalNumPrecond;i++)
	  zk1[I_precond[i]]+=V_precond[i]*rk[J_precond[i]];*/

	uint32_t rank;
	double tempV;
	/*for(uint32_t i = 0; i < dimension; ++i){
		zk1_compare[i] = 0;
		zk_compare[i] = 0;
	}
	
	for (uint32_t i = 0; i < dimension; ++i){
		tempV=0;
		a=row_idxL[i];
		b=row_idxL[i+1];
		for (uint32_t j=a;j<b;j++){
			tempV+=V_precond[j]*rk[J_precond[j]];
		}
		zk1_compare[i]+=tempV;
		//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
	}*/
	#pragma omp parallel private(rank,tempV, a, b)
	{
		rank=omp_get_thread_num();
		for (uint32_t i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
			tempV=0;
			a=row_idxL[i];
			b=row_idxL[i+1];
			for (uint32_t j=a;j<b;j++){
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
		for (uint32_t i=precondP_boundary[rank];i<precondP_boundary[rank+1];i++){
			tempV=0;
			a=row_idxLP[i];
			b=row_idxLP[i+1];
			for (uint32_t j=a;j<b;j++)
				tempV+=V_precondP[j]*zk1[J_precondP[j]];
			zk[i]+=tempV;
		}

	}				

	//initialize_all(dimension,pk_d,bp_d,x_d,zk_d,vector_in_d);
	#pragma omp parallel for
	for (uint32_t i=0;i<dimension;i++){
		pk[i]=zk[i];
		bp[i]=0;
		vector_out[i]=0;
	}


	//cublasDdot(handle, dimension,vector_in_d,1,vector_in_d,1,&doth);
	#pragma omp parallel for reduction(+:doth)
	for (uint32_t i=0;i<dimension;i++){
		doth+=vector_in[i]*vector_in[i];
	}
	double errorNorm=sqrt(doth);
	//printf("errorNorm is %f\n",errorNorm);
	uint32_t iter=0;
	
	struct timeval start1, end1;
	struct timeval start2, end2;

	gettimeofday(&start1, NULL);
	/*-----------Start the iteration by a while loop-----------*/
	while (iter<MAXIter&&error>1E-8){
		//matrix-vector multiplication, accelerated by our own functions
		/*#pragma omp parallel for
		  for (uint32_t i=0;i<totalNum_1;i++)
		  bp[I_1[i]]+=V_1[i]*pk[J_1[i]];*/
		#pragma omp parallel private(rank,tempV, a, b)
		{
			rank=omp_get_thread_num();
			//gettimeofday(&start2, NULL);
			for (uint32_t i=boundary[rank];i<boundary[rank+1];i++){
				tempV=0;
				a=row_idx[i];
				b=row_idx[i+1];
				for (uint32_t j=a;j<b;j++){
					tempV+=V[j]*pk[J[j]];
				}
				bp[i]+=tempV;
			}
			//gettimeofday(&end2, NULL);
			//printf("at iter %d, rank is %d, boundary span is %d, time is %d us\n", iter, rank, boundary[rank+1]-boundary[rank], ((end2.tv_sec * 1000000 + end2.tv_usec)-(start2.tv_sec * 1000000 + start2.tv_usec)));
		}

		dotp0=0;
		#pragma omp parallel for reduction(+:dotp0)
		for (uint32_t i=0;i<dimension;i++)
			dotp0+=bp[i]*pk[i];

		//r_k*z_k
		dotrz0=0;
		#pragma omp parallel for reduction(+:dotrz0)
		for (uint32_t i=0;i<dimension;i++)
			dotrz0+=rk[i]*zk[i];

		alphak=dotrz0/dotp0;
		alphak_1=-alphak;

		//update x_k to x_k+1, x=x+alphak*p;
		//update r_k to r_k+1, r=r-alphak*bp;
		#pragma omp parallel for
		for (uint32_t i=0;i<dimension;i++){
			vector_out[i]+=alphak*pk[i];
			rk[i]+=alphak_1*bp[i];
			zk1[i]=0;
			zk[i]=0;
		}

		//SpMV: zk=inv(m)*rk
		#pragma omp parallel private(rank,tempV, a, b)
		{
			rank=omp_get_thread_num();
			for (uint32_t i=precond_boundary[rank];i<precond_boundary[rank+1];i++){
				tempV=0;
				a=row_idxL[i];
				b=row_idxL[i+1];
				for (uint32_t j=a;j<b;j++){
					tempV+=V_precond[j]*rk[J_precond[j]];
				}
				zk1[i]+=tempV;
				//if (rank==4&&i<precond_boundary[4]+64) printf("tempV is %f\n",tempV);
			}
		}		

		#pragma omp parallel private(rank,tempV, a, b)
		{
			rank=omp_get_thread_num();
			for (uint32_t i=precondP_boundary[rank];i<precondP_boundary[rank+1];i++){
				tempV=0;
				a=row_idxLP[i];
				b=row_idxLP[i+1];
				for (uint32_t j=a;j<b;j++)
					tempV+=V_precondP[j]*zk1[J_precondP[j]];
				zk[i]+=tempV;
			}

		}			

		//r_k+1 * z_k+1
		dotrz1=0;
		#pragma omp parallel for reduction(+:dotrz1)
		for (uint32_t i=0;i<dimension;i++)
			dotrz1+=rk[i]*zk[i];
		betak=dotrz1/dotrz0;
		//printf("dotp0 is %f, dotrz0 is %f dotrz1 is %f, betak is %f and alphak is %f at iter %d\n", dotp0, dotrz0,dotrz1,betak, alphak, iter);
		//p=r+gamak*p;
		#pragma omp parallel for
		for (uint32_t i=0;i<dimension;i++){
			pk[i]=zk[i]+betak*pk[i];
			bp[i]=0;
		}
		dotr0=0;
		#pragma omp parallel for reduction (+:dotr0)
		for (uint32_t i=0;i<dimension;i++)
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


