#include "fspai.h"

void insertSort(int *J, float *V, int num, int *outJ, float *outV)
{
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


void fspaiCPU(S *SInput)
/*void fspai(int *I, int *J, float *V, int *I_precond, int *J_precond, float *V_precond, const int maxRowNum, const int *numInRow, 
			const int *row_idx, const int *numInRowPrecond, const int *row_idxPrecond, const float *diag, int colStart, int colEnd)*/
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
	int *row_idx=SInput->row_idx;
	int *numInRowPrecond=SInput->numInRowPrecond;
	int *row_idxPrecond=SInput->row_idxPrecond;
	float *diag=SInput->diag;
	int colStart=SInput->colStart; 
	int colEnd=SInput->colEnd;
	int id=SInput->id;
	
	float *sortedV= (float *)malloc(maxRowNum*sizeof(float));
	int *sortedJ=(int *)malloc(maxRowNum*sizeof(int));
	
	int start,num;
	//printf("fspai start at %d end at %d maxRowNum is %d\n", colStart, colEnd,maxRowNum);
	//Sort the matrix, using insert sort algorithm, will try quick sort algorithm
	for (int i=colStart;i<colEnd;i++)
	{
		
		start=row_idx[i];
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
	float *subMatrix=(float *) malloc(maxRowNum*maxRowNum*sizeof(float));
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
			index1=J[row_idx[i]+j];
			if(index1>i)
			{
				//first time touch the column index larger than i
				if (start==0) start=j;
				colIndex=0;
				for (int k=0; k<numInRow[index1]; k++)
				{
					index2=row_idx[index1]+k;
					if (J[index2]>i)
					{
						//travel at both row, the row i is currently in colIndex+j
						tempCol=J[row_idx[i]+colIndex+start];
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
								tempCol=J[row_idx[i]+colIndex+start];
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
				xk[subRowIndex]=V[row_idx[i]+j];
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
			index3=p+row_idxPrecond[i];
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
