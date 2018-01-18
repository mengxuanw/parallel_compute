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
	double start0, start1, stop,start;
	int size = (p - mypid -1)/numprocess + 3;
	double *A = (double*)malloc(sizeof(double)*m*n*size*q);
	double *X = (double*)malloc(sizeof(double)*m*n*size);
	double *B = (double*)malloc(sizeof(double)*m*n*size);
	int i,j,k,l;
	if(mypid == 0)
	{
		//read file
		double *A0 = (double*)malloc(sizeof(double)*m*n*p*q);
		double *X0 = (double*)malloc(sizeof(double)*m*n*p);
		double *B0 = (double*)malloc(sizeof(double)*m*n*p);
		start0 = MPI_Wtime();
		FILE *fileA;
		if((fileA = fopen("A.dat","rb")) == NULL)
		{
			printf("fileA open fail\n");
			return 0;
		}
		for(i = 0; i < m*n*p*q; i ++)
		{
			fread(&A0[i],sizeof(double),1,fileA);
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
			fread(&X0[i],sizeof(double),1,fileX);
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
			fread(&B0[i],sizeof(double),1,fileB);
		}
		fclose(fileB);
		stop = MPI_Wtime();
		printf("read file time:%0.16f\n",stop - start0);


		int *temp = (int*)malloc(sizeof(int)*numprocess);
		for(i = 0; i < numprocess; i ++)
			temp[i] = (p - i - 1)/numprocess + 1;
		int *count = (int*)malloc(sizeof(int)*numprocess * 2);
		count[0] = -1;
		count[numprocess] = temp[0];
		for(i = 1; i < numprocess; i ++)
		{
			count[i] = count[i - 1] + temp[i-1];
			count[i + numprocess] = count[i - 1 + numprocess] + temp[i];
		}

		//send data
		start1 = MPI_Wtime();
		for(int t = 1; t < numprocess -1; t ++)
		{
			int tempsize = (p - t -1)/numprocess + 3;
			double *tempA = (double*)malloc(sizeof(double)*m*n*tempsize*q);
			double *tempX = (double*)malloc(sizeof(double)*m*n*tempsize);
			double *tempB = (double*)malloc(sizeof(double)*m*n*tempsize);
			for(i = 0; i < m; i ++)
			{
				for(j = 0; j < n; j ++)
				{
					for(k = 0; k < tempsize; k ++)
					{
						tempX[i*n*tempsize + j*tempsize + k] = X0[i*n*p + j*p + k + count[t]];
						tempB[i*n*tempsize + j*tempsize + k] = B0[i*n*p + j*p + k + count[t]];
						for(l = 0; l < q; l ++)
						{
							tempA[i*n*tempsize*q + j*tempsize*q + k*q + l] = A0[i*n*p*q + j*p*q + (k+count[t])*q + l];
						}
					}
				}
			}
			MPI_Send(tempA, m*n*tempsize*q, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
			MPI_Send(tempX, m*n*tempsize, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
			MPI_Send(tempB, m*n*tempsize, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
			free(tempA);
			free(tempX);
			free(tempB);
		}

		//data for process numprocess-1
		int t = numprocess-1;
		int tempsize = (p - t -1)/numprocess + 3;
		double *tempA = (double*)malloc(sizeof(double)*m*n*tempsize*q);
		double *tempX = (double*)malloc(sizeof(double)*m*n*tempsize);
		double *tempB = (double*)malloc(sizeof(double)*m*n*tempsize);
		for(i = 0; i < m; i ++)
		{
			for(j = 0; j < n; j ++)
			{
				for(k = 0; k < tempsize-1; k ++)
				{
					tempX[i*n*tempsize + j*tempsize + k] = X0[i*n*p + j*p + k + count[t]];
					tempB[i*n*tempsize + j*tempsize + k] = B0[i*n*p + j*p + k + count[t]];
					for(l = 0; l < q; l ++)
					{
						tempA[i*n*tempsize*q + j*tempsize*q + k*q + l] = A0[i*n*p*q + j*p*q + (k+count[t])*q + l];
					}
				}
				tempX[i*n*tempsize + j*tempsize + k] = 0;//X0[i*n*tempsize + j*tempsize + k + count[t]];
				tempB[i*n*tempsize + j*tempsize + k] = 0;//B0[i*n*tempsize + j*tempsize + k + count[t]];
				for(l = 0; l < q; l ++)
				{
					tempA[i*n*tempsize*q + j*tempsize*q + k*q + l] = 0;//A0[i*n*tempsize*q + j*tempsize*q + (k+count[t])*q + l];
				}
			}
		}
		MPI_Send(tempA, m*n*tempsize*q, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
		MPI_Send(tempX, m*n*tempsize, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
		MPI_Send(tempB, m*n*tempsize, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
		free(tempA);
		free(tempX);
		free(tempB);

		//data for process 0
		t = 0;
		for(i = 0; i < m; i ++)
		{
			for(j = 0; j < n; j ++)
			{
				k = 0;
				X[i*n*size + j*size + k] = 0;//X0[i*n*tempsize + j*tempsize + k + count[t]];
				B[i*n*size + j*size + k] = 0;//B0[i*n*tempsize + j*tempsize + k + count[t]];
				for(l = 0; l < q; l ++)
				{
					A[i*n*size*q + j*size*q + k*q + l] = 0;//A0[i*n*tempsize*q + j*tempsize*q + (k+count[t])*q + l];
				}
				for(k = 1; k < size; k ++)
				{
					X[i*n*size + j*size + k] = X0[i*n*p + j*p + k + count[t]];
					//printf("i:%d,j:%d,k:%d,X[i*n*size + j*size + k]:%f\n",i,j,k,X[i*n*size + j*size + k]);
					B[i*n*size + j*size + k] = B0[i*n*p + j*p + k + count[t]];
					for(l = 0; l < q; l ++)
					{
						A[i*n*size*q + j*size*q + k*q + l] = A0[i*n*p*q + j*p*q + (k+count[t])*q + l];
					}
				}
			}
		}

		stop = MPI_Wtime();
		printf("send data time:%0.16f\n",stop - start1);
		start1 = MPI_Wtime();
	}
	else
	{
		MPI_Recv(A, m * n * size * q, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(X, m * n * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(B, m * n * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	//compute
	
	double sum_0 =0, sum_1 = 0, sum_2 = 0, sum_3 = 0;
	for(i = 0; i < m; i ++)
	{
		sum_1 = 0;
		for( j = 0; j < n; j ++)
		{
			sum_2 = 0;
			for(k = 1; k < size-1; k ++)
			{
				sum_3 = 0;
				int sym = (m/2 + i)%m;
				sum_3 += A[i*n*size*q + j*size*q + k*q + 0] * X[i*n*size + j*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 1] * X[(i==0?m-1:i-1)*n*size + j*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 2] * X[(i==m-1?0:i+1)*n*size + j*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 3] * X[(j==0?sym:i)*n*size + (j==0?j:j-1)*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 4] * X[(j==n-1?sym:i)*n*size + (j==n-1?j:j+1)*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 5] * X[(j==n-1?sym:(i==m-1?0:i+1))*n*size + (j==n-1?j:j+1)*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 6] * X[(j==0?sym:(i==m-1?0:i+1))*n*size + (j==0?j:j-1)*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 7] * X[(j==0?sym:(i==0?m-1:i-1))*n*size + (j==0?j:j-1)*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 8] * X[(j==n-1?sym:(i==0?m-1:i-1))*n*size + (j==n-1?j:j+1)*size + k];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 9] * X[i*n*size + j*size + k-1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 10] * X[(i==0?m-1:i-1)*n*size + j*size + k-1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 11] * X[(i==m-1?0:i+1)*n*size + j*size + k-1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 12] * X[(j==0?sym:i)*n*size + (j==0?j:j-1)*size + k-1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 13] * X[(j==n-1?sym:i)*n*size + (j==n-1?j:j+1)*size + k-1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 14] * X[i*n*size + j*size + k+1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 15] * X[(i==0?m-1:i-1)*n*size + j*size + k+1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 16] * X[(i==m-1?0:i+1)*n*size + j*size + k+1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 17] * X[(j==0?sym:i)*n*size + (j==0?j:j-1)*size + k+1];
				sum_3 += A[i*n*size*q + j*size*q + k*q + 18] * X[(j==n-1?sym:i)*n*size + (j==n-1?j:j+1)*size + k+1];

				sum_3 -= B[i*n*size + j*size + k];
				sum_2 += sum_3 * sum_3;	
			}
			sum_1 += sum_2;
		}
		sum_0 += sum_1;
	}
	if(mypid == 0)
	{
		stop = MPI_Wtime();
		printf("compute time:%0.16f\n",stop - start1);
		double sum = 0;
		for(i = 1; i < numprocess; i ++)
		{
			MPI_Recv(&sum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sum_0 += sum;
		}
		printf("result: %.20f\n", sum_0);
	}
	else
	{
		MPI_Send(&sum_0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	free(A);
	free(B);
	free(X);
	MPI_Finalize();
	return 0;
}
