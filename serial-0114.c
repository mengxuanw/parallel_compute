#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>

int main()
{
	int mypid, numprocess;
	int namelen;
	char   processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Init(NULL, NULL);//MPI Initialize
	MPI_Comm_rank(MPI_COMM_WORLD, &mypid);//获得当前进程号
	MPI_Comm_size(MPI_COMM_WORLD, &numprocess);//获得进程个数
	MPI_Get_processor_name(processor_name,&namelen);
	//printf("Process %d of %d is on %s\n",mypid, numprocess, processor_name);
	int m = 360, n = 180, p = 38, q = 19;
	double start0,start1,stop;
	int i,j,k,l;
	double *A = (double*)malloc(sizeof(double)*m*n*p*q);
	double *X = (double*)malloc(sizeof(double)*m*n*p);
	double *B = (double*)malloc(sizeof(double)*m*n*p);

	start0 = MPI_Wtime();
	FILE *fileA;
	if((fileA = fopen("A.dat","rb")) == NULL)
	{
		printf("fileA open fail\n");
		return 0;
	}
	for(i = 0; i < m*n*p*q; i ++)
	{
		fread(&A[i],sizeof(double),1,fileA);
	}
	fclose(fileA);

	FILE *fileX;
	if((fileX = fopen("x0.dat","rb")) == NULL)
	{
		printf("fileX open fail\n");
		return 0;
	}
	for(i = 0; i < m*n*p; i ++)
	{
		fread(&X[i],sizeof(double),1,fileX);
	}
	fclose(fileX);

	FILE *fileB;
	if((fileB = fopen("b.dat","rb")) == NULL)
	{
		printf("fileB open fail\n");
		return 0;
	}
	for(i = 0; i < m*n*p; i ++)
	{
		fread(&B[i],sizeof(double),1,fileB);
	}
	fclose(fileB);
	stop = MPI_Wtime();
	printf("read file time: %0.16f\n", stop - start0);

	start1 = MPI_Wtime();
	double sum_0 =0, sum_1 = 0, sum_2 = 0, sum_3 = 0;
	for(i = 0; i < m; i ++)
	{
		sum_1 = 0;
		for( j = 0; j < n; j ++)
		{
			sum_2 = 0;
			for(k = 0; k < p; k ++)
			{
				sum_3 = 0;
				int sym = (m/2 + i)%m;
				sum_3 += A[i*n*p*q + j*p*q + k*q + 0] * X[i*n*p + j*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 1] * X[(i==0?m-1:i-1)*n*p + j*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 2] * X[(i==m-1?0:i+1)*n*p + j*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 3] * X[(j==0?sym:i)*n*p + (j==0?j:j-1)*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 4] * X[(j==n-1?sym:i)*n*p + (j==n-1?j:j+1)*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 5] * X[(j==n-1?sym:(i==m-1?0:i+1))*n*p + (j==n-1?j:j+1)*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 6] * X[(j==0?sym:(i==m-1?0:i+1))*n*p + (j==0?j:j-1)*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 7] * X[(j==0?sym:(i==0?m-1:i-1))*n*p + (j==0?j:j-1)*p + k];
				sum_3 += A[i*n*p*q + j*p*q + k*q + 8] * X[(j==n-1?sym:(i==0?m-1:i-1))*n*p + (j==n-1?j:j+1)*p + k];
				if(k != 0)
				{
					sum_3 += A[i*n*p*q + j*p*q + k*q + 9] * X[i*n*p + j*p + k-1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 10] * X[(i==0?m-1:i-1)*n*p + j*p + k-1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 11] * X[(i==m-1?0:i+1)*n*p + j*p + k-1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 12] * X[(j==0?sym:i)*n*p + (j==0?j:j-1)*p + k-1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 13] * X[(j==n-1?sym:i)*n*p + (j==n-1?j:j+1)*p + k-1];
				}
				if(k != p - 1)
				{
					sum_3 += A[i*n*p*q + j*p*q + k*q + 14] * X[i*n*p + j*p + k+1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 15] * X[(i==0?m-1:i-1)*n*p + j*p + k+1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 16] * X[(i==m-1?0:i+1)*n*p + j*p + k+1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 17] * X[(j==0?sym:i)*n*p + (j==0?j:j-1)*p + k+1];
					sum_3 += A[i*n*p*q + j*p*q + k*q + 18] * X[(j==n-1?sym:i)*n*p + (j==n-1?j:j+1)*p + k+1];
				}
				sum_3 -= B[i*n*p + j*p + k];
				sum_2 += sum_3 * sum_3;	
			}
			sum_1 += sum_2;
		}
		sum_0 += sum_1;
	}
	stop = MPI_Wtime();
	printf("compute time: %0.16f\n",stop - start1);
	printf("result: %.20f\n", sum_0);
	free(A);
	free(B);
	free(X);
	MPI_Finalize();
	return 0;
}