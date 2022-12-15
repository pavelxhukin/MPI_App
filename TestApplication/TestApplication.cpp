/*
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <ctime>

double f(double a) {
	return(4.0 / (1.0 + a * a));
}

int main(int argc, char** argv) {

	int done = 0, n, myid, numprocs, i;
	double PI25 = 3.141592653589793238462643;
	double mypi, pi, h, sum, x;
	double startwtime = 0.0, endwtime;
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Get_processor_name(processor_name, &namelen);

	fprintf(stderr, "Process %d on %s\n", myid, processor_name);

	n = 0;
	while (!done) {
		if (myid == 0) {
			//print("enter the number of intervals: (0 quits) ");
			//fflush(stdout);
			//scanf("%d",&n);

			if (n == 0) n = 100; else n = 0;

			startwtime = MPI_Wtime();
		}

		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (n == 0)
			done = 1;
		else {
			h = 1.0 / (double)n;
			sum = 0.0;
			for (i = myid + 1; i <= n; i += numprocs) {
				x = h * ((double)i - 0.5);
				sum += f(x);
			}
			mypi = h * sum;

			MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			if (myid == 0) {
				printf("pi is approximately %.16f, Error is %.16f\n",
					pi, fabs(pi - PI25));
				endwtime = MPI_Wtime();
				printf("wall clock time = %f\n",
					endwtime - startwtime);
			}
		}
	}

	MPI_Finalize();
	return 0;
}*/

// 1-22 ******************************************************************************************************************

/*
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstring>
#include <string>
#include "mpi.h"
using namespace std;

void DataInitialization(double* x, int N)
{
	for (int i = 0; i < N; i++)
	{
		*(x + i) = (rand() % 110) - 99;
	}
}

void DataPrint(double* x, int N)
{
	for (int i = 0; i < N; i++)
	{
		cout << "[" << i << "]" << *(x + i) << endl;
	}
}

int main(int argc, char* argv[]) {
	double* x = 0;
	int	TotalMulti = 0;
	int ProcRank, ProcNum, N = 0, K = 0;
	MPI_Status Status;
	// инициализация
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Barrier(MPI_COMM_WORLD);
	double time = MPI_Wtime();
	// подготовка данных
	if (ProcRank == 0) {
		cout << "Write K: "; 
		cin >> K;
		cout << "N:";
		cin >> N;
	}
	MPI_Bcast(&K, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	x = new double[N];
	if (ProcRank == 0)
	{
		DataInitialization(x, N);
		DataPrint(x, N);
	}
	// рассылка данных на все процессы

	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//обработка данных на процессах
	int k = N / ProcNum;
	int count = 0;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) { i2 = N; }
	for (int i = i1; i < i2; i++)
	{
		if (x[i] > 0)
		{
			if (x[i]==K)
			count++;
		}
	}
	// сборка частичных сумм на процессе с рангом 0 
	MPI_Reduce(&count,&TotalMulti, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	// вывод результата 
	if (ProcRank == 0)
	{
		cout << " \nKol: " << TotalMulti << endl;
		double timeend = MPI_Wtime();
		double past = timeend - time;
		cout << "Time: " << past << " (seconds?)" << endl;
	}
	delete[] x;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
} */


//1 - 23 ******************************************************************************************************************


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstring>
#include <string>
#include "mpi.h"
using namespace std;

void DataInitialization(int* x, int N)
{
	for (int i = 0; i < N; i++)
	{
		*(x + i) = rand() % 300;
	}
}

void DataPrint(int* x, int N)
{
	for (int i = 0; i < N; i++)
	{
		cout << "[" << i << "]" << *(x + i) << endl;
	}
}

int main(int argc, char* argv[]) {
	int* x = 0;
	int	TotalMulti = 0;
	int ProcRank, ProcNum, N = 0;
	MPI_Status Status;
	// инициализация
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Barrier(MPI_COMM_WORLD);
	double time = MPI_Wtime();
	// подготовка данных
	if (ProcRank == 0) {
		cout << "N:";
		cin >> N;
	}
	MPI_Bcast(&N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	x = new int[N];
	if (ProcRank == 0)
	{
		DataInitialization(x, N);
		DataPrint(x, N);
	}
	// рассылка данных на все процессы

	MPI_Bcast(x, N, MPI_INT, 0, MPI_COMM_WORLD);

	//обработка данных на процессах
	int k = N / ProcNum;
	int sum = 0;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) { i2 = N; }
	for (int i = i1; i < i2; i++)
	{
		if (!((i+1) % 3))
		{
			sum+= x[i];
		}
	}
	// сборка частичных сумм на процессе с рангом 0 
	MPI_Reduce(&sum, &TotalMulti, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	// вывод результата 
	if (ProcRank == 0)
	{
		cout << " \nKol: " << TotalMulti << endl;
		double timeend = MPI_Wtime();
		double past = timeend - time;
		cout << "Time: " << past << " (seconds?)" << endl;
	}
	delete[] x;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

