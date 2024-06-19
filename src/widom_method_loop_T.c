//gcc -o widom_method_loop_T widom_method_loop_T.c -lgsl -lgslcblas -lm
//./widom_method_loop_T 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

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
                e += 4 * (r6i * r6i - r6i) - (shift ? ecut : 0.0);
                *vir += 48 * (r6i * r6i - 0.5 * r6i);
            }
        }
    }
    return e + (tailcorr ? ecor : 0.0);
}

double total_e(double *rx, double *ry, double *rz, int N, double L,
               double rc2, int tailcorr, double ecor,
               int shift, double ecut, double *vir) {
    int i;
    double tvir;
    double e = 0.0;
    *vir = 0.0;
    for (i = 0; i < N - 1; i++) {
        e += energy_inter(i, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &tvir, i + 1);
        *vir += tvir;
    }
    return e;
}

void write_xyz(FILE *fp, double *rx, double *ry, double *rz, int n, double L) {
    int i;
    fprintf(fp, "%i\n", n);
    fprintf(fp, "BOX %.5lf %.5lf %.5lf\n", L, L, L);
    for (i = 0; i < n; i++) {
        fprintf(fp, "%.5lf %.5lf %.5lf\n", rx[i], ry[i], rz[i]);
    }
}

void init(double *rx, double *ry, double *rz, int n, double L, gsl_rng *r) {
    int i, ix, iy, iz;
    int n3 = 2;
    while ((n3 * n3 * n3) < n) n3++;
    ix = iy = iz = 0;
    for (i = 0; i < n; i++) {
        rx[i] = ((double)ix + 0.5) * L / n3;
        ry[i] = ((double)iy + 0.5) * L / n3;
        rz[i] = ((double)iz + 0.5) * L / n3;
        ix++;
        if (ix == n3) {
            ix = 0;
            iy++;
            if (iy == n3) {
                iy = 0;
                iz++;
            }
        }
    }
}

void widom(double *rx, double *ry, double *rz, int N, double L,
           double rc2, int shift, double ecut, gsl_rng *r, double *e) {
    int j;
    double dx, dy, dz, hL = L / 2.0, r2, r6i;
    (*e) = 0.0;
    rx[N] = (gsl_rng_uniform(r) - 0.5) * L;
    ry[N] = (gsl_rng_uniform(r) - 0.5) * L;
    rz[N] = (gsl_rng_uniform(r) - 0.5) * L;
    for (j = 0; j < N; j++) {
        dx = (rx[N] - rx[j]);
        dy = (ry[N] - ry[j]);
        dz = (rz[N] - rz[j]);
        if (dx > hL) dx -= L;
        else if (dx < -hL) dx += L;
        if (dy > hL) dy -= L;
        else if (dy < -hL) dy += L;
        if (dz > hL) dz -= L;
        else if (dz < -hL) dz += L;
        r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < rc2) {
            r6i = 1.0 / (r2 * r2 * r2);
            *e += 4 * (r6i * r6i - r6i) - (shift ? ecut : 0.0);
        }
    }
}

enum { XYZ, NONE };

int main(int argc, char *argv[]) {
    double *rx, *ry, *rz;
    int N = 216, c, a, p;
    double L = 0.0;
    double we, w_sum;
    double beta;
    double rho = 0.5, T_start = 100.0, T_end = 700, T_step = 50, T, rc2 = 3.5, vir, vir_old, p_sum, pcor, V;
    double E_new, E_old, esum, rr3, ecor, ecut;
    double ei_new, ei_old, ivir_new, ivir_old;
    double dr = 0.2, dx, dy, dz;
    double rxold, ryold, rzold;
    int i, j;
    int nCycles = 10, nSamp, nEq = 1000;
    int nAcc;
    int short_out = 0;
    int shift = 0;
    int tailcorr = 1;
    int prog = 0;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long int Seed = 23410981;
    FILE *fp;
    char *traj_fn = NULL;
    int traj_out = XYZ;
    int traj_samp = 100;
    double sum_e_minus_beta_W = 0.0;
    int nInsertions = 0;

    // Open CSV file for writing
    FILE *csv_fp = fopen("temperature_output.csv", "w");
    if (csv_fp == NULL) {
        perror("Unable to open temperature_output.csv for writing");
        return 1;
    }

    // Write the CSV header
    fprintf(csv_fp, "T,Average Energy per Particle,Average Pressure,Excess Chemical Potential\n");

    for (T = T_start; T <= T_end; T += T_step) {
        rc2 = 3.5;
        gsl_rng_set(r, Seed);
        L = pow(((double)N / rho), 1.0 / 3.0);
        V = L * L * L;
        beta = 1.0 / T;
        ecor = (tailcorr ? (8.0 * M_PI * rho / 3.0) * (1.0 / 3.0 * pow(rc2, -3.0) - pow(rc2, -1.5)) : 0.0);
        ecut = (4.0 * (pow(rc2, -6.0) - pow(rc2, -3.0)));
        nSamp = nCycles * N;

        rx = (double *)malloc(N * sizeof(double));
        ry = (double *)malloc(N * sizeof(double));
        rz = (double *)malloc(N * sizeof(double));
        init(rx, ry, rz, N, L, r);
        esum = 0.0;
        p_sum = 0.0;
        nAcc = 0;
        fp = (traj_fn && traj_out != NONE) ? fopen(traj_fn, "w") : NULL;

        for (i = -nEq; i < nSamp; i++) {
            j = (int)(gsl_rng_uniform(r) * N);
            rxold = rx[j];
            ryold = ry[j];
            rzold = rz[j];
            E_old = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &vir_old, 0);
            rx[j] += (gsl_rng_uniform(r) - 0.5) * dr;
            ry[j] += (gsl_rng_uniform(r) - 0.5) * dr;
            rz[j] += (gsl_rng_uniform(r) - 0.5) * dr;
            if (rx[j] < -0.5 * L) rx[j] += L;
            else if (rx[j] > 0.5 * L) rx[j] -= L;
            if (ry[j] < -0.5 * L) ry[j] += L;
            else if (ry[j] > 0.5 * L) ry[j] -= L;
            if (rz[j] < -0.5 * L) rz[j] += L;
            else if (rz[j] > 0.5 * L) rz[j] -= L;
            E_new = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &vir, 0);
            if (gsl_rng_uniform(r) < exp(-beta * (E_new - E_old))) {
                nAcc++;
                vir_old = vir;
            } else {
                rx[j] = rxold;
                ry[j] = ryold;
                rz[j] = rzold;
                E_new = E_old;
            }

            if (i >= 0) {
                esum += E_new;
                p_sum += vir_old / (3.0 * V);
                if (traj_fn && traj_out != NONE && (i % traj_samp == 0)) {
                    write_xyz(fp, rx, ry, rz, N, L);
                }
                widom(rx, ry, rz, N, L, rc2, shift, ecut, r, &we);
                sum_e_minus_beta_W += exp(-beta * we);
                nInsertions++;
            }
        }

        esum /= (double)nSamp;
        p_sum = rho * T + p_sum / (double)nSamp;
        if (traj_fn && traj_out != NONE) fclose(fp);

        // Compute excess chemical potential
        double excess_chemical_potential = -T * log(sum_e_minus_beta_W / nInsertions);

        // Write results to CSV
        fprintf(csv_fp, "%.2f,%.5f,%.5f,%.5f\n", T, esum / N, p_sum, excess_chemical_potential);

        free(rx);
        free(ry);
        free(rz);
    }

    fclose(csv_fp);
    gsl_rng_free(r);

    return 0;
}
