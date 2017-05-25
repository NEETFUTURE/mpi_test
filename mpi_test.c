#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define N 100
#define NSTEP 600
#define UMU 1.257e-6
#define EPS0 8.854e-12
#define C 2.998e8
#define SIGMA 0.0
#define PI 3.141592
#define freq 0.5e9
#define STEP 1

int main(int argc, char **argv)
{
    int n, k, myid, numprocs, i;
    double recvbuf, sendbuf;

    double h[N];
    double e[N];
    double dt, ec1, ec2, hc;
    double t;
    double dz = 1.0e-2;

    char filename[20];

    MPI_File ffile;

    MPI_Request req1, req2;
    MPI_Status stat1, stat2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    //printf("Program Start!\n");

    t = 0.0;
    dt = dz / C;
    ec1 = (1.0 - SIGMA * dt / (2.0 * EPS0)) / (1.0 + SIGMA * dt / (2.0 * EPS0));
    ec2 = -(dt / (EPS0 * dz) / (1.0 + SIGMA * dt / (2.0 * EPS0)));
    hc = -dt / (dz * UMU);

    //printf("In rank%d,  ec1 = %f, ec2 = %f, hc = %f\n", myid,  ec1, ec2, hc);

    for (i = 0; i < N; i++)
    {
        e[i] = 0.0;
        h[i] = 0.0;
    }


    for (n = 0; n < NSTEP; n++)
    {


        // ファイルの保存
        sprintf(filename, "data_mpi_test/data%05d.raw", n);

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
            sendbuf = e[1];
            MPI_Isend(&sendbuf, 1, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, &req1);
            MPI_Irecv(&recvbuf, 1, MPI_DOUBLE, myid - 1, 1, MPI_COMM_WORLD, &req1);
            MPI_Waitall(1,&req1, &stat1);
            e[0] = recvbuf;
        }
        if (myid < numprocs - 1)
        {
            sendbuf = e[N - 2];
            MPI_Isend(&sendbuf, 1, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, &req2);
            MPI_Irecv(&recvbuf, 1, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, &req2);
            MPI_Waitall(1,&req2, &stat2);
            e[N - 1] = recvbuf;
        }

        printf("LOOPING\n");
        MPI_Barrier(MPI_COMM_WORLD);

        if (myid == 1 && t < 0.5 / freq)
        {
            e[N / 2] = e[N / 2] + pow(sin(2.0 * PI * freq * t), 4);
        }

        for (k = 1; k < N - 1; k++)
        {
            e[k] = ec1 * e[k] + ec2 * (h[k] - h[k - 1]);
        }

        t = t + dt / 2.0;

        for (k = 0; k < N - 1; k++)
        {
            h[k] = h[k] + hc * (e[k + 1] - e[k]);
        }

        t = t + dt / 2.0;
    }

    MPI_File_close(&ffile);

    MPI_Finalize();
}
