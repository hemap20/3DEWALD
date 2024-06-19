//gcc -o widom_method widom_method.c -lgsl -lgslcblas -lm
//./widom_method -N 216 -rho 0.5 -T 1.0 -dr 0.2 -rc 3.5 -nc 10 -ne 1000 +tc -sh -prog 100 -s 23410981 -traj_samp 100 -traj output.xyz
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
//sum of all pair interactions between a particle i and all particles from i0 to N-1.
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
//writes the particle configuration to a file in XYZ format
void write_xyz(FILE *fp, double *rx, double *ry, double *rz, int n, double L) {
    int i;
    fprintf(fp, "%i\n", n);
    fprintf(fp, "BOX %.5lf %.5lf %.5lf\n", L, L, L);
    for (i = 0; i < n; i++) {
        fprintf(fp, "%.5lf %.5lf %.5lf\n", rx[i], ry[i], rz[i]);
    }
}
//initializes the particle positions by placing them on a cubic grid and then scaling them to fit within the simulation box of size L.
void init(double *rx, double *ry, double *rz, int n, double L, gsl_rng *r) {
    int i, ix, iy, iz;
    int n3 = 2;
    //Finds the smallest cube of grid points that can fit all particles
    while ((n3 * n3 * n3) < n) n3++;
    ix = iy = iz = 0;
    for (i = 0; i < n; i++) {
      //Arrays of the x, y, and z coordinates of particles
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
        //if the distance is within rc2, the Lennard-Jones potential is calculated and added to the total energy *e
        if (r2 < rc2) {
            r6i = 1.0 / (r2 * r2 * r2);
            *e += 4 * (r6i * r6i - r6i) - (shift ? ecut : 0.0);
        }
    }
}

enum { XYZ, NONE };
// initializes variables, sets up the random number generator, reads input parameters
// initializes particle positions, and performs Monte Carlo simulations
int main(int argc, char *argv[]) {
    double *rx, *ry, *rz;
    int N = 216, c, a, p;
    double L = 0.0;
    double we, w_sum;
    double beta;
    double rho = 0.5, T = 1.0, rc2 = 3.5, vir, vir_old, p_sum, pcor, V;
    double E_new, E_old, esum, rr3, ecor, ecut;
    double ei_new, ei_old, ivir_new, ivir_old;
    double dr = 0.2, dx, dy, dz; //dr = max displacement for MC simulations
    double rxold, ryold, rzold;
    int i, j;
    int nCycles = 10, nSamp, nEq = 1000; //nSamp = The total number of sampling steps, calculated as nCycles * N.
    int nAcc; //number of accepted moves
    int short_out = 0;
    int shift = 0; //A flag for applying a shift in the potential energy calculation, not detailed in the provided segment.
    int tailcorr = 1; //A flag indicating whether to apply the tail correction to the potential energy.
    int prog = 0;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long int Seed = 23410981;
    FILE *fp;
    char *traj_fn = NULL;
    int traj_out = XYZ;
    int traj_samp = 100; //The frequency of trajectory sampling for output.
    double sum_e_minus_beta_W = 0.0;
    int nInsertions = 0; //The number of particle insertions in the Widom test for chemical potential calculation.

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N")) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rho")) rho = atof(argv[++i]);
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
    /* For the first nEq steps, the system equilibrates. After that, the system 
    samples configurations and calculates the energy and pressure. The loop includes 
    attempts to move a random particle and applies the Metropolis criterion to decide 
    whether to accept the move. It also handles periodic boundary conditions.*/
    for (i = -nEq; i < nSamp; i++) {
        // Select a random particle to move
        j = (int)(gsl_rng_uniform(r) * N);
        rxold = rx[j];
        ryold = ry[j];
        rzold = rz[j];
        // Calculate the old energy of the selected particle
        E_old = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &vir_old, 0);
        // Attempt to move the particle
        rx[j] += (gsl_rng_uniform(r) - 0.5) * dr;
        ry[j] += (gsl_rng_uniform(r) - 0.5) * dr;
        rz[j] += (gsl_rng_uniform(r) - 0.5) * dr;
        // Apply periodic boundary conditions
        if (rx[j] < -0.5 * L) rx[j] += L;
        else if (rx[j] > 0.5 * L) rx[j] -= L;
        if (ry[j] < -0.5 * L) ry[j] += L;
        else if (ry[j] > 0.5 * L) ry[j] -= L;
        if (rz[j] < -0.5 * L) rz[j] += L;
        else if (rz[j] > 0.5 * L) rz[j] -= L;
        // Calculate the new energy of the selected particle
        E_new = energy_inter(j, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &vir, 0);
        // Metropolis criterion to decide whether to accept the move
        if (gsl_rng_uniform(r) < exp(beta * (E_old - E_new))) {
            nAcc++;
            vir_old = vir;
        } else {
            // Revert the move if not accepted
            rx[j] = rxold;
            ry[j] = ryold;
            rz[j] = rzold;
        }
        // Sampling phase
        if (i >= 0) {
            esum += E_new;
            p_sum += vir_old;
            if ((i % N) == 0) {
                // Output progress and trajectory data
                if (prog > 0 && ((i / N) % prog) == 0) {
                    printf("Cycle %i of %i\n", i / N - nEq, nCycles);
                }
                if (traj_fn && traj_out == XYZ && ((i / N) % traj_samp) == 0) {
                    write_xyz(fp, rx, ry, rz, N, L);
                }
                // Widom insertion method to calculate excess chemical potential
                widom(rx, ry, rz, N, L, rc2, shift, ecut, r, &we);
                sum_e_minus_beta_W += exp(-beta * we);
                nInsertions++;
            }
        }
    }
    esum /= nSamp;
    p_sum = rho * T + p_sum / (3.0 * V * nSamp);
    if (traj_fn) fclose(fp);

    // Calculate the excess chemical potential
    double average_e_minus_beta_W = sum_e_minus_beta_W / nInsertions;
    double mu_ex = -T * log(average_e_minus_beta_W / rho);

    printf("E/NkT = %.5lf\n", esum / (N * T));
    printf("P/rho kT = %.5lf\n", p_sum / rho);
    printf("Excess chemical potential (mu_ex) = %.5lf\n", mu_ex);
    printf("Number of particle insertions = %.5d\n", nInsertions);

    free(rx);
    free(ry);
    free(rz);
    gsl_rng_free(r);
    return 0;
}
