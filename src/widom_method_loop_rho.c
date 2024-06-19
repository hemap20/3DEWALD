//gcc -o widom_method_loop_rho widom_method_loop_rho.c -lgsl -lgslcblas -lm
//./widom_method_loop_rho -N 216 -rho_start 0.3 -rho_end 0.8 -rho_step 0.1 -T 1.0 -dr 0.2 -rc 3.5 -nc 10000 -ne 1000
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
    double rho_start = 0.3, rho_end = 0.8, rho_step = 0.1; // Define the range and step for rho
    double T = 1.0, rc2 = 3.5, vir, vir_old, p_sum, pcor, V;
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

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N")) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rho_start")) rho_start = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rho_end")) rho_end = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rho_step")) rho_step = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T")) T = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dr")) dr = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rc")) rc2 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-nc")) nCycles = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-ne")) nEq = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s")) Seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-traj_samp")) traj_samp = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-traj")) {
            traj_fn = argv[++i];
            traj_out = XYZ;
        }
        else if (!strcmp(argv[i], "-prog")) prog = atoi(argv[++i]);
    }

    rc2 *= rc2;
    gsl_rng_set(r, Seed);
    beta = 1.0 / T;

    rx = (double *)malloc(N * sizeof(double));
    ry = (double *)malloc(N * sizeof(double));
    rz = (double *)malloc(N * sizeof(double));

    for (double rho = rho_start; rho <= rho_end; rho += rho_step) {
        L = pow(((double)N / rho), 1.0 / 3.0);
        V = L * L * L;
        ecor = (tailcorr ? (8.0 * M_PI * rho / 3.0) * (1.0 / 3.0 * pow(rc2, -3.0) - pow(rc2, -1.5)) : 0.0);
        ecut = (4.0 * (pow(rc2, -6.0) - pow(rc2, -3.0)));
        pcor = (tailcorr ? (32.0 * M_PI * rho * rho / 9.0) * (pow(rc2, -3.0) - 1.5 * pow(rc2, -1.5)) : 0.0);
        
        init(rx, ry, rz, N, L, r);
        esum = 0.0;
        w_sum = 0.0;
        p_sum = 0.0;
        nAcc = 0;
        nSamp = 0;
        fp = NULL;
        if (traj_fn != NULL) fp = fopen(traj_fn, "w");

        E_old = total_e(rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &vir);
        for (i = 0; i < nEq; i++) {
            for (j = 0; j < N; j++) {
                rxold = rx[j];
                ryold = ry[j];
                rzold = rz[j];
                ei_old = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &ivir_old, 0);
                rx[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                ry[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                rz[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                if (rx[j] > L) rx[j] -= L;
                if (rx[j] < 0) rx[j] += L;
                if (ry[j] > L) ry[j] -= L;
                if (ry[j] < 0) ry[j] += L;
                if (rz[j] > L) rz[j] -= L;
                if (rz[j] < 0) rz[j] += L;
                ei_new = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &ivir_new, 0);
                if (gsl_rng_uniform(r) < exp(-beta * (ei_new - ei_old))) {
                    E_old += (ei_new - ei_old);
                    vir += (ivir_new - ivir_old);
                    nAcc++;
                } else {
                    rx[j] = rxold;
                    ry[j] = ryold;
                    rz[j] = rzold;
                }
            }
            if (prog) printf("Equilibration step %d/%d\r", i + 1, nEq);
        }

        for (i = 0; i < nCycles; i++) {
            for (j = 0; j < N; j++) {
                rxold = rx[j];
                ryold = ry[j];
                rzold = rz[j];
                ei_old = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &ivir_old, 0);
                rx[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                ry[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                rz[j] += (gsl_rng_uniform(r) - 0.5) * dr;
                if (rx[j] > L) rx[j] -= L;
                if (rx[j] < 0) rx[j] += L;
                if (ry[j] > L) ry[j] -= L;
                if (ry[j] < 0) ry[j] += L;
                if (rz[j] > L) rz[j] -= L;
                if (rz[j] < 0) rz[j] += L;
                ei_new = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &ivir_new, 0);
                if (gsl_rng_uniform(r) < exp(-beta * (ei_new - ei_old))) {
                    E_old += (ei_new - ei_old);
                    vir += (ivir_new - ivir_old);
                    nAcc++;
                } else {
                    rx[j] = rxold;
                    ry[j] = ryold;
                    rz[j] = rzold;
                }
            }
            E_new = E_old;
            esum += E_new;
            p_sum += rho * T + (vir / (3.0 * V)) + pcor;

            for (p = 0; p < traj_samp; p++) {
                widom(rx, ry, rz, N, L, rc2, shift, ecut, r, &we);
                sum_e_minus_beta_W += exp(-beta * we);
                nInsertions++;
            }
            if (fp && (i % traj_samp) == 0) {
                write_xyz(fp, rx, ry, rz, N, L);
            }
            nSamp++;
            if (prog) printf("Cycle %d/%d\r", i + 1, nCycles);
        }

        if (fp) fclose(fp);

        double avgE = esum / (double)nSamp;
        double avgP = p_sum / (double)nSamp;
        double mu_excess = -log(sum_e_minus_beta_W / (double)nInsertions) / beta;

        printf("rho = %.2f, Average Energy per Particle = %.5lf, Average Pressure = %.5lf, Excess Chemical Potential = %.5lf\n", rho, avgE / N, avgP, mu_excess);
    }

    free(rx);
    free(ry);
    free(rz);
    gsl_rng_free(r);
    return 0;
}
