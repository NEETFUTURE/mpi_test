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
    double recvbuf1, sendbuf1, recvbuf2, sendbuf2;

    double h[N];
    double e[N];
    double dt, ec1, ec2, hc;
    double t;
    double dz = 1.0e-2;


    MPI_Request req1, req2;
    MPI_Status stat1, stat2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    t = 0.0;
    dt = dz / C;
    ec1 = (1.0 - SIGMA * dt / (2.0 * EPS0)) / (1.0 + SIGMA * dt / (2.0 * EPS0));
    ec2 = -(dt / (EPS0 * dz) / (1.0 + SIGMA * dt / (2.0 * EPS0)));
    hc = -dt / (dz * UMU);


    for (i = 0; i < N; i++)
    {
        e[i] = 0.0;
        h[i] = 0.0;
    }

    for (n = 0; n < NSTEP; n++)
    {




        if (myid == 1)
        {
            sendbuf1 = h[0];
            printf("id=%d, n=%d, send = %8.5f\n", myid, n, sendbuf1);
            fflush(stdout);
            MPI_Isend(&sendbuf1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &req1);
            MPI_Irecv(&recvbuf1, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &req1);
            MPI_Waitall(1,&req1, &stat1);
            e[0] = recvbuf1;
        }
        if (myid == 0)
        {
            sendbuf2 = e[N - 1];
            MPI_Isend(&sendbuf2, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &req2);
            MPI_Irecv(&recvbuf2, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &req2);
            MPI_Waitall(1,&req2, &stat2);
            printf("id=%d, n=%d, recv = %8.5f\n", myid, n, recvbuf2);
            printf("======================================\n");
            fflush(stdout);

        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (myid == 1 && t < 0.5 / freq)
        {
            e[N / 2] = e[N / 2] + pow(sin(2.0 * PI * freq * t), 4);
        }

        for (k = 1; k < N; k++)
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

    MPI_Finalize();
}
