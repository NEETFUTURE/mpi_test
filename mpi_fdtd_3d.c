#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define N 30
#define NSTEP 80
#define UMU 1.257e-6
#define EPS0 8.854e-12
#define C 2.998e8
#define SIGMA 0.0
#define PI 3.141592
#define freq 0.6e9
#define STEP 1

void write_to_file(double *e, double *h, FILE *fp_p, int myid)
{
    int i;

    if (myid == 0)
    {
        for (i = N - 4; i < N; i++)
        {
            fprintf(fp_p, "%9.7f        ", e[i]);
        }
        fprintf(fp_p, "\n");
        for (i = N - 4; i < N; i++)
        {
            fprintf(fp_p, "        %9.7f", h[i]);
        }
        fprintf(fp_p, "\n");
    }
    else
    {
        for (i = 0; i < 4; i++)
        {
            fprintf(fp_p, "%9.7f        ", e[i]);
        }
        fprintf(fp_p, "\n");
        for (i = 0; i < 4; i++)
        {
            fprintf(fp_p, "        %9.7f", h[i]);
        }
        fprintf(fp_p, "\n");
    }
}

int main(int argc, char **argv)
{
    int n, k, myid, numprocs, i, j;
    double recvbuf1, sendbuf1, recvbuf2, sendbuf2;

    double hx[N];
    double hy[N];
    double hz[N];

    double ex[N];
    double ey[N];
    double ez[N];

    double dt, ec1, ec2, hc;
    double t;
    double dz = 1.0e-2;

    char filename[20];
    char filename_p[20];

    MPI_File ffile;
    FILE *fp_p;

    MPI_Request req1, req2;
    MPI_Status stat1, stat2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    t = 0.0;
    dt = dz / C;
    ec1 = (1.0 - SIGMA * dt / (2.0 * EPS0)) / (1.0 + SIGMA * dt / (2.0 * EPS0));
    ec2 = (dt / (EPS0 * dz) / (1.0 + SIGMA * dt / (2.0 * EPS0)));
    hc = -dt / (dz * UMU);


    for (i = 0; i < N; i++)
    {
        ex[i] = ey[i] = ez[i] = 0.0;
        hx[i] = hy[i] = hz[i] = 0.0;
    }


    sprintf(filename_p, "data_mpi_test_p/data_%d.txt", myid);
    fp_p = fopen(filename_p, "w");

    for (n = 0; n < NSTEP; n++)
    {
        // ファイルの保存
        sprintf(filename, "data_mpi_3d/data%05d.raw", n);

        write_to_file(e, h, fp_p, myid);
        fprintf(fp_p, "----------------------------------------------------------------------\n",n);

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

        // ******************* 電界の計算 *********************

        if (myid == 0 && t < 0.5 / freq)
        {
            ez[N/2] = ez[N/2] + pow(sin(2.0 * PI * freq * t), 4);
        }

        for (i = 0; i < N; i++)

        {
            for (j = 0; j < N; j++)
            {
                for (k = 1; k < N; k++)
                {
                    ex[i][j][k] = ec1 * ex[i][j][k] + ec2 * (-hy[i][j][k] + hy[i][j][k - 1] + hy[i][j][k] - hy[i][j-1][k]);
                    ey[i][j][k] = ec1 * ey[i][j][k] + ec2 * ( hx[i][j][k] - hx[i][j][k - 1] - hz[i][j][k] + hz[i-1][j][k]);
                    ez[i][j][k] = ec1 * ez[i][j][k] + ec2 * (-hx[i][j][k] + hx[i][j-1][k] + hy[i][j][k] - hy[i-1][j][k]);
                }
            }
        }

        write_to_file(e, h, fp_p, myid);
        fprintf(fp_p, "----------------------------------------------------------------------\n",n);

        // ******************* 電界の送受信 *********************
        if (myid == 0){
            sendbuf2 = e[N - 1];
            MPI_Isend(&sendbuf2, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &req1);
        } else {
            MPI_Irecv(&recvbuf1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &req2);
            MPI_Waitall(1,&req2, &stat2);
            e[0] = recvbuf1;
        }

        t = t + dt / 2.0;

        write_to_file(e, h, fp_p, myid);
        fprintf(fp_p, "----------------------------------------------------------------------\n",n);

        // ******************* 磁界の計算 *********************
        for (k = 0; k < N - 1; k++)
        {
            h[k] = h[k] + hc * (e[k + 1] - e[k]);
        }

        t = t + dt / 2.0;

        write_to_file(e, h, fp_p, myid);
        fprintf(fp_p, "----------------------------------------------------------------------\n",n);

        // ******************* 磁界の送受信 *********************
        if (myid == 1){
            sendbuf1 = h[0];
            MPI_Isend(&sendbuf1, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &req2);
        } else {
            MPI_Irecv(&recvbuf2, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &req1);
            MPI_Waitall(1,&req1, &stat1);
            h[N - 1] = recvbuf2;
        }

        write_to_file(e, h, fp_p, myid);

        MPI_Barrier(MPI_COMM_WORLD);


        fprintf(fp_p, "n=%04d**************************************************************\n",n);
    }

    MPI_File_close(&ffile);
    fclose(fp_p);

    MPI_Finalize();
}
