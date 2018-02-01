#include "fspai.h"

/*float for preconditioner generating only*/
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

void solverPrecondCPU(const int dimension, const int totalNum, const int *I_accum, const int *J, const float *V, const int totalNumPrecond, const int *I_precond_accum,
			const int *J_precond, const float *V_precond, const int totalNumPrecondP, const int *I_precondP_accum, const int *J_precondP, const float *V_precondP, 
			 const float *y, float *x, const int MAXIter, int *realIter)
{
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
???LINES MISSING
???LINES MISSING
???LINES MISSING
???LINES MISSING
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
