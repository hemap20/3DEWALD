//gcc -o widom_method_loop_rho widom_method_loop_rho.c -lgsl -lgslcblas -lm
//./widom_method_loop_rho
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

// Function prototypes
double energy_inter(int i, double *rx, double *ry, double *rz, int N, double L,
                    double rc2, int tailcorr, double ecor,
                    int shift, double ecut, double *vir, int i0);

double total_e(double *rx, double *ry, double *rz, int N, double L,
               double rc2, int tailcorr, double ecor,
               int shift, double ecut, double *vir);

void init(double *rx, double *ry, double *rz, int n, double L, gsl_rng *r);

void widom(double *rx, double *ry, double *rz, int N, double L,
           double rc2, int shift, double ecut, gsl_rng *r, double *e);

void write_xyz(FILE *fp, double *rx, double *ry, double *rz, int n, double L);

enum { XYZ, NONE };

int main(int argc, char *argv[]) {
    double *rx, *ry, *rz;
    int N = 216;
    double L = 0.0;
    double beta;
    double rho_start = 0.3, rho_end = 0.8, rho_step = 0.1;
    double T = 1.0, rc2 = 3.5, vir, pcor, V;
    double E_new, E_old, esum, ecor, ecut;
    double ei_new, ei_old, ivir_new, ivir_old;
    double dr = 0.2;
    double rxold, ryold, rzold;
    int i, j;
    int nCycles = 10, nSamp, nEq = 1000;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long int Seed = 23410981;
    FILE *fp;
    double sum_e_minus_beta_W;
    int nInsertions;
    int traj_samp = 100;

    // Open CSV file for writing results
    fp = fopen("output.csv", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening output.csv for writing.\n");
        return 1;
    }

    // Write headers to CSV file
    fprintf(fp, "rho,Average Energy per Particle,Average Pressure,Excess Chemical Potential\n");

    for (double rho = rho_start; rho <= rho_end; rho += rho_step) {
        L = pow(((double)N / rho), 1.0 / 3.0);
        V = L * L * L;
        ecor = 8.0 * M_PI * rho / 3.0 * (1.0 / 3.0 * pow(rc2, -3.0) - pow(rc2, -1.5));
        ecut = 4.0 * (pow(rc2, -6.0) - pow(rc2, -3.0));
        pcor = 32.0 * M_PI * rho * rho / 9.0 * (pow(rc2, -3.0) - 1.5 * pow(rc2, -1.5));

        rx = (double *)malloc(N * sizeof(double));
        ry = (double *)malloc(N * sizeof(double));
        rz = (double *)malloc(N * sizeof(double));

        init(rx, ry, rz, N, L, r);
        esum = 0.0;
        sum_e_minus_beta_W = 0.0;
        nInsertions = 0;
        nSamp = 0;
        E_old = total_e(rx, ry, rz, N, L, rc2, 1, ecor, 0, ecut, &vir);

        // Equilibration phase
        for (i = 0; i < nEq; i++) {
            for (j = 0; j < N; j++) {
                rxold = rx[j];
                ryold = ry[j];
                rzold = rz[j];
                ei_old = energy_inter(j, rx, ry, rz, N, L, rc2, 1, ecor, 0, ecut, &ivir_old, 0);
                rx[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                ry[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                rz[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                if (rx[j] > L) rx[j] -= L;
                if (rx[j] < 0) rx[j] += L;
                if (ry[j] > L) ry[j] -= L;
                if (ry[j] < 0) ry[j] += L;
                if (rz[j] > L) rz[j] -= L;
                if (rz[j] < 0) rz[j] += L;
                ei_new = energy_inter(j, rx, ry, rz, N, L, rc2, 1, ecor, 0, ecut, &ivir_new, 0);
                if (gsl_rng_uniform(r) < exp(-beta * (ei_new - ei_old))) {
                    E_old += (ei_new - ei_old);
                    vir += (ivir_new - ivir_old);
                } else {
                    rx[j] = rxold;
                    ry[j] = ryold;
                    rz[j] = rzold;
                }
            }
        }

        // Production phase
        for (i = 0; i < nCycles; i++) {
            for (j = 0; j < N; j++) {
                rxold = rx[j];
                ryold = ry[j];
                rzold = rz[j];
                ei_old = energy_inter(j, rx, ry, rz, N, L, rc2, 1, ecor, 0, ecut, &ivir_old, 0);
                rx[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                ry[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                rz[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                if (rx[j] > L) rx[j] -= L;
                if (rx[j] < 0) rx[j] += L;
                if (ry[j] > L) ry[j] -= L;
                if (ry[j] < 0) ry[j] += L;
                if (rz[j] > L) rz[j] -= L;
                if (rz[j] < 0) rz[j] += L;
                ei_new = energy_inter(j, rx, ry, rz, N, L, rc2, 1, ecor, 0, ecut, &ivir_new, 0);
                if (gsl_rng_uniform(r) < exp(-beta * (ei_new - ei_old))) {
                    E_old += (ei_new - ei_old);
                    vir += (ivir_new - ivir_old);
                    sum_e_minus_beta_W += exp(-beta * (ei_new - E_old));
                    nInsertions++;
                } else {
                    rx[j] = rxold;
                    ry[j] = ryold;
                    rz[j] = rzold;
                }
            }
            nSamp++;
        }

        double avgE = E_old / (double)(nCycles * N);
        double avgP = rho * T + vir / (3.0 * V) + pcor;
        double mu_excess = -log(sum_e_minus_beta_W / (double)nInsertions) / beta;

        // Write results to CSV file
        fprintf(fp, "%.5lf,%.5lf,%.5lf,%.5lf\n", rho, avgE, avgP, mu_excess);

        free(rx);
        free(ry);
        free(rz);
    }

    gsl_rng_free(r);
    fclose(fp);

    return 0;
}

// Function definitions

double energy_inter(int i, double *rx, double *ry, double *rz, int N, double L,
                    double rc2, int tailcorr, double ecor,
                    int shift, double ecut, double *vir, int i0) {
    int j;
    double dx, dy, dz, r2, r6i;
    double e = 0.0, hL = L / 2.0;
    *vir = 0.0;
    for (j = i0; j < N; j++) {
        if (i != j) {
            dx = rx[i] - rx[j];
            dy = ry[i] - ry[j];
            dz = rz[i] - rz[j];
            if (dx > hL) dx -= L;
            else if (dx < -hL) dx += L;
            if (dy > hL) dy -= L;
            else if (dy < -hL) dy += L;
            if (dz > hL) dz -= L;
            else if (dz < -hL) dz += L;
            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc2) {
                r6i = 1.0 / (r2 * r2 * r2);
                e += 4.0 * (r6i * r6i - r6i) - ecut;
                *vir += 24.0 * r6i * (2.0 * r6i - 1.0);
            }
        }
    }
    return e + ecor;
}

double total_e(double *rx, double *ry, double *rz, int N, double L,
               double rc2, int tailcorr, double ecor,
               int shift, double ecut, double *vir) {
    int i, j;
    double e = 0.0, vir_tmp;
    *vir = 0.0;
    for (i = 0; i < N; i++) {
        for (j = i + 1; j < N; j++) {
            double dx = rx[i] - rx[j];
            double dy = ry[i] - ry[j];
            double dz = rz[i] - rz[j];
            if (dx > L / 2.0) dx -= L;
            else if (dx < -L / 2.0) dx += L;
            if (dy > L / 2.0) dy -= L;
            else if (dy < -L / 2.0) dy += L;
            if (dz > L / 2.0) dz -= L;
            else if (dz < -L / 2.0) dz += L;
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc2) {
                double r6i = 1.0 / (r2 * r2 * r2);
                e += 4.0 * (r6i * r6i - r6i) - ecut;
                *vir += 24.0 * r6i * (2.0 * r6i - 1.0);
            }
        }
    }
    return e * 0.5 + ecor;
}

void init(double *rx, double *ry, double *rz, int n, double L, gsl_rng *r) {
    int i;
    double s3 = pow(n * 1.0, 1.0 / 3.0), hL = L / 2.0;
    int nx = (int)(s3), ny = (int)(s3), nz = (int)(s3);
    double dx = L / (nx + 1), dy = L / (ny + 1), dz = L / (nz + 1);
    int ntot = nx * ny * nz;
    if (ntot < n) {
        fprintf(stderr, "Density too low\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        rx[i] = dx * (1 + i % nx) + gsl_rng_uniform(r);
        ry[i] = dy * (1 + (i / nx) % ny) + gsl_rng_uniform(r);
        rz[i] = dz * (1 + i / (nx * ny)) + gsl_rng_uniform(r);
        if (rx[i] > L) rx[i] -= L;
        if (ry[i] > L) ry[i] -= L;
        if (rz[i] > L) rz[i] -= L;
    }
}

void write_xyz(FILE *fp, double *rx, double *ry, double *rz, int n, double L) {
    int i;
    fprintf(fp, "%d\n", n);
    fprintf(fp, "L = %lf\n", L);
    for (i = 0; i < n; i++) {
        fprintf(fp, "C  %15.7lf %15.7lf %15.7lf\n", rx[i], ry[i], rz[i]);
    }
}
