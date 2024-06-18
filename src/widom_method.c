#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

double energy_inter (int i, double * rx, double * ry, double * rz, int N, double L,
                     double rc2, int tailcorr, double ecor,
                     int shift, double ecut, double * vir, int i0) {
  int j;
  double dx, dy, dz, r2, r6i;
  double e = 0.0, hL = L / 2.0;
  *vir = 0.0;
  for (j = i0; j < N; j++) {
    if (i != j) {
      dx  = rx[i] - rx[j];
      dy  = ry[i] - ry[j];
      dz  = rz[i] - rz[j];
      if (dx > hL) dx -= L;
      else if (dx < -hL) dx += L;
      if (dy > hL) dy -= L;
      else if (dy < -hL) dy += L;
      if (dz > hL) dz -= L;
      else if (dz < -hL) dz += L;
      r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < rc2) {
        r6i   = 1.0 / (r2 * r2 * r2);
        e    += 4 * (r6i * r6i - r6i) - (shift ? ecut : 0.0);
        *vir += 48 * (r6i * r6i - 0.5 * r6i);
      }
    }
  }
  return e + (tailcorr ? ecor : 0.0);
}

double total_e (double * rx, double * ry, double * rz, int N, double L,
                double rc2, int tailcorr, double ecor,
                int shift, double ecut, double * vir) {
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

void write_xyz(FILE * fp, double * rx, double * ry, double * rz, int n, double L) {
  int i;
  fprintf(fp, "%i\n", n);
  fprintf(fp, "BOX %.5lf %.5lf %.5lf\n", L, L, L);
  for (i = 0; i < n; i++) {
    fprintf(fp, "%.5lf %.5lf %.5lf\n", rx[i], ry[i], rz[i]);
  }
}

void init (double * rx, double * ry, double * rz, int n, double L, gsl_rng * r) {
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

void widom (double * rx, double * ry, double * rz, int N, double L,
            double rc2, int shift, double ecut, gsl_rng * r, double * e) {
  int j;
  double dx, dy, dz, hL = L / 2.0, r2, r6i;
  (*e) = 0.0;
  rx[N] = (gsl_rng_uniform(r) - 0.5) * L;
  ry[N] = (gsl_rng_uniform(r) - 0.5) * L;
  rz[N] = (gsl_rng_uniform(r) - 0.5) * L;
  for (j = 0; j < N; j++) {
    dx  = (rx[N] - rx[j]);
    dy  = (ry[N] - ry[j]);
    dz  = (rz[N] - rz[j]);
    if (dx > hL) dx -= L;
    else if (dx < -hL) dx += L;
    if (dy > hL) dy -= L;
    else if (dy < -hL) dy += L;
    if (dz > hL) dz -= L;
    else if (dz < -hL) dz += L;
    r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < rc2) {
      r6i   = 1.0 / (r2 * r2 * r2);
      *e    += 4 * (r6i * r6i - r6i) - (shift ? ecut : 0.0);
    }
  }
}

enum { XYZ, NONE };

int main (int argc, char * argv[]) {
  double * rx, * ry, * rz;
  int N = 216, c, a, p;
  double L = 0.0;
  double we, w_sum;
  double beta;
  double rho = 0.5, T = 1.0, rc2 = 3.5, vir, vir_old, p_sum, pcor, V;
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
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;
  FILE * fp;
  char * traj_fn = NULL;
  int traj_out = XYZ;
  int traj_samp = 100;
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
  pcor = (tailcorr ? (16.0 * M_PI * rho * rho / 3.0) * (2.0 / 3.0 * pow(rc2, -3.0) - pow(rc2, -1.5)) : 0.0);
  ecut = (shift ? 4 * (pow(rc2, -6.0) - pow(rc2, -3.0)) : 0.0);
  rx = (double *)malloc((N + 1) * sizeof(double));
  ry = (double *)malloc((N + 1) * sizeof(double));
  rz = (double *)malloc((N + 1) * sizeof(double));
  init(rx, ry, rz, N, L, r);
  esum = 0.0;
  p_sum = 0.0;
  nSamp = 0;
  if (traj_out == XYZ) {
    if (traj_fn == NULL) fp = stdout;
    else {
      fp = fopen(traj_fn, "w");
      if (fp == NULL) {
        fprintf(stderr, "Could not open %s for writing\n", traj_fn);
        return 1;
      }
    }
  }
  for (c = -nEq; c < nCycles; c++) {
    for (a = 0; a < N; a++) {
      i = gsl_rng_uniform_int(r, N);
      rxold = rx[i];
      ryold = ry[i];
      rzold = rz[i];
      ei_old = energy_inter(i, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &ivir_old, 0);
      dx = (gsl_rng_uniform(r) - 0.5) * dr;
      dy = (gsl_rng_uniform(r) - 0.5) * dr;
      dz = (gsl_rng_uniform(r) - 0.5) * dr;
      rx[i] += dx;
      ry[i] += dy;
      rz[i] += dz;
      if (rx[i] < 0.0) rx[i] += L;
      else if (rx[i] >= L) rx[i] -= L;
      if (ry[i] < 0.0) ry[i] += L;
      else if (ry[i] >= L) ry[i] -= L;
      if (rz[i] < 0.0) rz[i] += L;
      else if (rz[i] >= L) rz[i] -= L;
      ei_new = energy_inter(i, rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &ivir_new, 0);
      if (exp(beta * (ei_old - ei_new)) > gsl_rng_uniform(r)) {
        vir += ivir_new - ivir_old;
        esum += ei_new - ei_old;
      } else {
        rx[i] = rxold;
        ry[i] = ryold;
        rz[i] = rzold;
      }
    }
    if (c >= 0) {
      esum += total_e(rx, ry, rz, N, L, rc2, tailcorr, ecor, shift, ecut, &vir);
      w_sum = 0.0;
      for (j = 0; j < N; j++) {
        widom(rx, ry, rz, N, L, rc2, shift, ecut, r, &we);
        w_sum += exp(-beta * we);
      }
      nSamp++;
      p_sum += vir / (3.0 * V) + rho * T;
      if (traj_out == XYZ && (c % traj_samp == 0)) write_xyz(fp, rx, ry, rz, N, L);
    }
    if (prog && ((c + nEq) % prog == 0)) fprintf(stderr, "Cycle %i of %i\n", c, nCycles);
  }
  if (traj_out == XYZ && traj_fn != NULL) fclose(fp);
  printf("E/NkT = %.5lf\n", esum / nSamp / N / T);
  printf("P/rho kT = %.5lf\n", p_sum / nSamp / rho / T);
  gsl_rng_free(r);
  free(rx);
  free(ry);
  free(rz);
  return 0;
}
