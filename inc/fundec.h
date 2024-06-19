///declare all functions here


double duration(struct timeb start, struct timeb end);

void print_lammps_input_file(double **PosIons, float *chg, int natoms, float **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep, char printmode, string filename);

void print_coor(double **PosIons, int natoms, float **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep,char printmode, string filename);

void print_carcoor(double **PosIons, int natoms, float **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep,char printmode, string filename);

void printCoor(double **PosIons, int natoms, string *type);

void printFor(double **ForceIons, int natoms, string *type);

void printVel(double **Vel, int natoms, string *type);

void printMD(double **PosIons, int natoms, float **boxcell, double **Vel, int MDstep,int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **fixatoms, int **pairs, float *mass);

void pairlist(double **PosIons, int natoms, float **boxcell, string *atomtype, int *natoms_type, int n_atomtype, string *type);


void readbondpar(int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int **pairs, int natoms);

void readnonbondpar(int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **pairs, int natoms);

double CalEner(double **PosIons, int natoms, float **boxcell, int nbatoms, int *atom1, int *atom2, string *npot, double **npar);

double CalEner1(double **PosIons, int natoms, float **boxcell, int *batom1, int *batom2, string *bondpot, int *nbatom1, int *nbatom2, string *nonbondpot, double **bondpar, double **nonbondpar, int **pairs, int randatom, int intdoub, int randatom2);

void CalFor(double **PosIons, double **ForIons,  int natoms, float **boxcell, int nbatoms, int *atom1, int *atom2, string *npot, double **npar);

double CalEnerFor(double **PosIons, int natoms, float **boxcell, double **ForceIons, string filename, bool bonded, bool forces);

double CalKinEner(double **Vel, int natoms, float *mass);

float CalTemp(double KEner, int natoms);

void res_vel(double **vel, int natoms, float Temp, float *mass);


void wrap_coor(double **PosIons, int natoms, float **boxcell);

void montecarlo(double **PosIons, int natoms, float **boxcell, float *mass, string *type, int n_atomtype, int *natoms_type, string *atomtype,float Temp, int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **fixatoms, int **pairs);


void minimization(double **PosIons, int natoms, double **ForceIons, float **boxcell, double Pot, string *type, int **fixatoms);


pair <double, float>  first_wolfe_brac(double **PosIons, double **PrevPos,double **ForceIons, int natoms,  float **boxcell, double Pot, float c1_wolfe, float alpha, float brac_par, int **fixatoms);


void md(double **PosIons, int natoms, float **boxcell, double **vel, double **ForceIons, float *mass, float dt, float Temp, int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **fixatoms);


void simuAnn(double **PosIons, int natoms, float **boxcell, double **vel, double **ForceIons, float *mass, float dt, float TempB, float TempE, int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **fixatoms, int rampsize, int rampstep);

void CalDens(double **PosIons, int natoms, float **boxcell, int n_atomtype, int *natoms_type, double lim1, double dx, double *density);

double selfe(int n_atomtype, int *natoms_type, float *chargs, float betaa);

double compute_real_energy(double **atomPositions, float *atomCharges, int numAtoms, double beta, float **simulationBox)

double dist(double **PosIons, int atom1, int atom2, float **box);

double reci_energy_h(int K, float **box, float *ion_charge, int n_atoms, float **pos_ions, double beta)

// void scalarProductMat(double mat[3][3], double k);

long double M_n(long double u, int n);

complex<long double> B(int m, int n, int K);

long double dotProduct(long double *v1, long double *v2);

long double dotProductu(double *v1, long double *v2);

void crossProduct(float *v_A, float *v_B, long double *out);

double bspline(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int M, int n);

long double error(long double a, long double b);