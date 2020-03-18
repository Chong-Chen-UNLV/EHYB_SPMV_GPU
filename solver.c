#include "kernel.h"
#include "solver.h"
#include "convert.h"



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

