#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define N 50
#define NSTEP 600
#define STEP 1

int main(int argc, char **argv)
{
    int n, k, myid, numprocs, i;
    double recvbuf1, sendbuf1, recvbuf2, sendbuf2;

    double e[N];

    char filename[20];

    MPI_File ffile;

    MPI_Request req1, req2;
    MPI_Status stat1, stat2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    for (i = 0; i < N; i++)
    {
        e[i] = 0.0;
	if (i == N/2 && myid == 0)
	    e[i] = 1.0;
    }


    for (n = 0; n < NSTEP; n++)
    {


        // ファイルの保存
        sprintf(filename, "data_mpi_test_test/data%05d.raw", n);

        MPI_File_open(
            MPI_COMM_WORLD, filename,
            MPI_MODE_WRONLY | MPI_MODE_CREATE,
            MPI_INFO_NULL, &ffile
        );

        MPI_File_set_view(
            ffile, myid*(N-1)*sizeof(double),
            MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL
        );
        MPI_File_write(ffile, e, N-1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        // ファイルの保存 END

        if (myid > 0)
        {
            sendbuf1 = e[1];
            MPI_Isend(&sendbuf1, 1, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, &req1);
            MPI_Irecv(&recvbuf1, 1, MPI_DOUBLE, myid - 1, 1, MPI_COMM_WORLD, &req1);
            MPI_Waitall(1,&req1, &stat1);
            e[0] = recvbuf1;
        }
        if (myid < numprocs - 1)
        {
            sendbuf2 = e[N - 2];
            MPI_Isend(&sendbuf2, 1, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, &req2);
            MPI_Irecv(&recvbuf2, 1, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, &req2);
            MPI_Waitall(1,&req2, &stat2);
            e[N - 1] = recvbuf2;
        }

        MPI_Barrier(MPI_COMM_WORLD);


        for (k = 1; k < N-1; k++)
        {
            e[k] = abs(e[k-1] - e[k+1]);
        }

    }

    MPI_File_close(&ffile);

    MPI_Finalize();
}
