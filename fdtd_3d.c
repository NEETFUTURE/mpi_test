#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 21
#define NY 21
#define NZ 21
#define NSTEP 300
#define UMU 1.257e-6
#define EPS0 8.854e-12
#define C 2.998e8
#define SIGMA 0.0
#define PI 3.141592
#define freq 0.6e9
#define STEP 1

#define TRGT_X 10
#define TRGT_Y 10
#define TRGT_Z 10

int main(int argc, char **argv)
{
    int n, i, j, k;

    double *hx = malloc(sizeof(double) * NX * NY * NZ);
    double *hy = malloc(sizeof(double) * NX * NY * NZ);
    double *hz = malloc(sizeof(double) * NX * NY * NZ);

    double *ex = malloc(sizeof(double) * NX * NY * NZ);
    double *ey = malloc(sizeof(double) * NX * NY * NZ);
    double *ez = malloc(sizeof(double) * NX * NY * NZ);

    double dt, ec1, ec2, hc;
    double t;
    double dz = 1.0e-2;

    char filename[20];

    FILE *sample_fp;
    FILE *fp;

    t = 0.0;
    dt = dz / C/ 3.0;
    ec1 = (1.0 - SIGMA * dt / (2.0 * EPS0)) / (1.0 + SIGMA * dt / (2.0 * EPS0));
    ec2 = dt / (EPS0 * dz) / (1.0 + SIGMA * dt / (2.0 * EPS0));
    hc = -dt / (dz * UMU);

    for (i = 0; i < NX; i++)
    {
        for (j = 0; j < NY; j++)
        {
            for (k = 0; k < NZ; k++)
            {
                ex[k*NX*NY + NX*j + i] = 0.0;
                ey[k*NX*NY + NX*j + i] = 0.0;
                ez[k*NX*NY + NX*j + i] = 0.0;
                hx[k*NX*NY + NX*j + i] = 0.0;
                hy[k*NX*NY + NX*j + i] = 0.0;
                hz[k*NX*NY + NX*j + i] = 0.0;
            }
        }
    }


    // sample_fp = fopen("sample.csv", "w");
    //fprintf(sample_fp, "time,e\n");

    for (n = 0; n < NSTEP; n++)
    {

        sprintf(filename, "data_3d/data_%04d.csv", n);
        fp = fopen(filename, "w");
        fprintf(fp, "x,y,z,e\n");


        // ******************* 電界の計算 *********************
        if (t < 0.5 / freq)
        {
            ex[TRGT_X*NX*NY + TRGT_Y*NX + TRGT_Z] = ex[TRGT_X*NX*NY + TRGT_Y*NX + TRGT_Z] + pow(sin(2.0 * PI * freq * (t+dz*i/C)), 4);
            printf("%9.7f\n", pow(sin(2.0 * PI * freq * (t+dz*i/C)), 4));
        }

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                for (k = 0; k < NZ; k++)
                {
                    if (j > 0 && k > 0)
                    {
                        ex[k*NX*NY + NX*j + i] = ec1 * ex[k*NX*NY + NX*j + i]
                                               + ec2 * (- hy[k*NX*NY + NX*j + i]
                                                        + hy[k*NX*NY + NX*j + i - 1]
                                                        + hz[k*NX*NY + NX*j + i]
                                                        - hz[k*NX*NY + NX*(j - 1) + i]);
                    }
                    if (i > 0 && k > 0)
                    {
                        ey[k*NX*NY + NX*j + i] = ec1 * ey[k*NX*NY + NX*j + i]
                                               + ec2 * (+ hx[k*NX*NY + NX*j + i]
                                                        - hx[k*NX*NY + NX*j + i - 1]
                                                        - hz[k*NX*NY + NX*j + i]
                                                        + hz[(k-1)*NX*NY + NX*j + i]);
                    }
                    if (i > 0 && j > 0){
                        ez[k*NX*NY + NX*j + i] = ec1 * ez[k*NX*NY + NX*j + i]
                                               + ec2 * (- hx[k*NX*NY + NX*j + i]
                                                        + hx[k*NX*NY + NX*(j-1) + i]
                                                        + hy[k*NX*NY + NX*j + i]
                                                        - hy[(k-1)*NX*NY + NX*j + i]);
                    }
                    fprintf(fp, "%d, %d, %d, %9.7f\n", i, j, k, ex[k*NX*NY + NX*j + i]);

                }
            }
        }

        //fprintf(sample_fp, "%d,%9.7f\n", n, ex[25*N*N + N*14 + 14]);

        t = t + dt / 2.0;

        // ******************* 磁界の計算 *********************
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                for (k = 0; k < NZ; k++)
                {
                    if (j < NY-1 && k < NZ-1){
                        hx[k*NX*NY + NX*j + i] = hx[k*NX*NY + NX*j + i]
                                               + hc * (+ ey[k*NX*NY + NX*j + i]
                                                       - ey[k*NX*NY + NX*j + i + 1]
                                                       - ez[k*NX*NY + NX*j + i]
                                                       + ez[k*NX*NY + NX*(j+1) + i]);
                    }

                    if (i < NX-1 && k < NZ-1){
                        hy[k*NX*NY + NX*j + i] = hy[k*NX*NY + NX*j + i]
                                               + hc * (- ex[k*NX*NY + NX*j + i]
                                                       + ex[k*NX*NY + NX*j + i + 1]
                                                       + ez[k*NX*NY + NX*j + i]
                                                       - ez[(k+1)*NX*NY + NX*j + i]);
                    }

                    if (j < NY-1 && i < NX-1){
                        hz[k*NX*NY + NX*j + i] = hz[k*NX*NY + NX*j + i]
                                               + hc * (+ ex[k*NX*NY + NX*j + i]
                                                       - ex[k*NX*NY + NX*(j+1) + i]
                                                       - ey[k*NX*NY + NX*j + i]
                                                       + ey[(k+1)*NX*NY + NX*j + i]);
                    }
                }
            }
        }

        t = t + dt / 2.0;
        fclose(fp);
    }
    // fclose(sample_fp);


}
