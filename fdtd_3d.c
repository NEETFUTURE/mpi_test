#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 200
#define NY 200
#define NZ 5
#define NSTEP 300
#define UMU 1.257e-6
#define EPS0 8.854e-12
#define C 2.998e8
#define SIGMA 0.0
#define PI 3.141592
#define freq 0.6e9
#define STEP 1

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

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
            {
                ex[i*NX*NY + NX*j + k] = 0.0;
                ey[i*NX*NY + NX*j + k] = 0.0;
                ez[i*NX*NY + NX*j + k] = 0.0;
                hx[i*NX*NY + NX*j + k] = 0.0;
                hy[i*NX*NY + NX*j + k] = 0.0;
                hz[i*NX*NY + NX*j + k] = 0.0;
            }
        }
    }


    sample_fp = fopen("sample.csv", "w");
    fprintf(sample_fp, "time,e\n");

    for (n = 0; n < NSTEP; n++)
    {
        printf("n = %04d\n", n);

        sprintf(filename, "data_3d/data_%04d.csv", n);
        fp = fopen(filename, "w");

        fprintf(fp, "x,y,z,e\n");


        // ******************* 電界の計算 *********************
        if (t < 0.5 / freq)
        {
            ex[40*NX*NY + NX*(NY/2) + (NX/2)] = ex[i*NX*NY + NX*(NY/2) + (NX/2)] + pow(sin(2.0 * PI * freq * (t+dz*i/C)), 4);
            ex[40*NX*NY + NX*(NY/2) + (NX/2)] = ex[i*NX*NY + NX*(NY/2) + (NX/2)] + pow(sin(2.0 * PI * freq * (t+dz*i/C)), 4);
            ey[40*NX*NY + NX*(NY/2) + (NX/2)] = ex[i*NX*NY + NX*(NY/2) + (NX/2)] + pow(sin(2.0 * PI * freq * (t+dz*i/C)), 4);
            ey[41*NX*NY + NX*(NY/2) + (NX/2)] = ex[i*NX*NY + NX*(NY/2) + (NX/2)] + pow(sin(2.0 * PI * freq * (t+dz*i/C)), 4);
        }

        for (i = 1; i < N; i++)
        {
            for (j = 1; j < N; j++)
            {
                for (k = 1; k < N; k++)
                {
                    if (i < NX - 1)
                    {
                        ex[i*NX*NY + NX*j + k] = ec1 * ex[i*NX*NY + NX*j + k]
                                               + ec2 * (- hy[i*NX*NY + NX*j + k]
                                                        + hy[i*NX*NY + NX*j + k - 1]
                                                        + hz[i*NX*NY + NX*j + k]
                                                        - hz[i*NX*NY + NX*(j - 1) + k]);
                    }
                    if (j < NY - 1)
                    {
                        ey[i*NX*NY + NX*j + k] = ec1 * ey[i*NX*NY + NX*j + k]
                                               + ec2 * (+ hx[i*NX*NY + NX*j + k]
                                                        - hx[i*NX*NY + NX*j + k - 1]
                                                        - hz[i*NX*NY + NX*j + k]
                                                        + hz[(i-1)*NX*NY + NX*j + k]);
                    }
                    if (k < NZ - 1){
                        ez[i*NX*NY + NX*j + k] = ec1 * ez[i*NX*NY + NX*j + k]
                                               + ec2 * (- hx[i*NX*NY + NX*j + k]
                                                        + hx[i*NX*NY + NX*(j-1) + k]
                                                        + hy[i*NX*NY + NX*j + k]
                                                        - hy[(i-1)*NX*NY + NX*j + k]);
                    }
                    fprintf(fp, "%d, %d, %d, %9.7f\n", i, j, k, ex[i*NX*NY + NX*j + k]);
                }
            }
        }
        fprintf(sample_fp, "%d,%9.7f\n", n, ex[25*N*N + N*25 + 25]);

        t = t + dt / 2.0;

        // ******************* 磁界の計算 *********************
        for (i = 0; i < N - 1; i++)
        {
            for (j = 0; j < N - 1; j++)
            {
                for (k = 0; k < N - 1; k++)
                {
                    hx[i*NX*NY + NX*j + k] = hx[i*NX*NY + NX*j + k]
                                           + hc * (+ ey[i*NX*NY + NX*j + k]
                                                   - ey[i*NX*NY + NX*j + k + 1]
                                                   - ez[i*NX*NY + NX*j + k]
                                                   + ez[i*NX*NY + NX*(j+1) + k]);

                    hy[i*NX*NY + NX*j + k] = hy[i*NX*NY + NX*j + k]
                                           + hc * (- ex[i*NX*NY + NX*j + k]
                                                   + ex[i*NX*NY + NX*j + k + 1]
                                                   + ez[i*NX*NY + NX*j + k]
                                                   - ez[(i+1)*NX*NY + NX*j + k]);

                    hz[i*NX*NY + NX*j + k] = hz[i*NX*NY + NX*j + k]
                                           + hc * (+ ex[i*NX*NY + NX*j + k]
                                                   - ex[i*NX*NY + NX*(j+1) + k]
                                                   - ey[i*NX*NY + NX*j + k]
                                                   + ey[(i+1)*NX*NY + NX*j + k]);
                }
            }
        }

        t = t + dt / 2.0;
        fclose(fp);
    }
    fclose(sample_fp);


}
