//calculation of reciprocal energy
#include "libinclude.h"
#include "const.h"
#include <vector>
//initialising reciprocal space vector



//4 implementations for the energy calculation

#define BASIC 1
#define REDUCTION_REAL_IMG 2
#define REDUCTION_KVECTOR 3

//creating reciprocal space vectors
#if defined BASIC
    double reci_energy(int K, float **box, float *ion_charge, int n_atoms, float **pos_ions, double beta){
        double reci_energy_h= 0
        for(int kx=-K; kx<K+1; kx++){
            for (int ky=-K; ky<K+1; ky++){
                for(int kz=-K; kz<K+1; kz++){
                    //assigining the vectors
                    std::vector<double> G(3);
                    G[0] = 2*M_PI*kx/box[0][0];
                    G[1] = 2*M_PI*ky/box[1][1];
                    G[2] = 2*M_PI*kz/box[2][2];
                    complex<double> t(0,1);
                    complex<double> SG = 0;
                    for(int i=0, i<n_atoms, i++){
                        //SG is calculated for every possible reciprocal vector summed over the number of atoms for better accuracy  
                        double G_dot_r = G[0]*pos_ion[i][0] + G[1]*pos_ion[i][1] + G[2]*pos_ion[i][2];
                        SG += ion_charge[i]*(cos(G_dot_r) + t*sin(G_dot_r));
                    }
                    double norm_SG = norm(SG);
                    double mod_G = sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]);
                    double damping_factor_exp = -(mod_G*mod_G)/(4*beta*beta);
                    reci_energy_h+= (exp(damping_factor_exp)*norm_SG)/(mod_G*mod_G);
                }
            }
        }
        //normalisation of energy
        reci_energy_h= reci_energy*2*M_PI/(box[0][0]*box[1][1]*box[2][2]);
        return reci_energy;
    }
#elif defined REDUCTION_REAL_IMG
//separate loops for real and imaginary
    double reci_energy(int K, float **box, float *ion_charge, int n_atoms, float **pos_ions, double beta){
        double reci_energy_h= 0;
        complex<double> t(0,1);
        complex<double> SG = 0;
        //iterating through the reciprocal vectors
        for(int kx<-K; kx<K+1; kx++){
            for(int ky=-K; ky<K+1; ky++){
                for(int kz=-K; kz<K+1; kz++){
                    std::vector<double> G(3);
                    G[0] = 2*M_PI*kx/box[0][0];
                    G[1] = 2*M_PI*ky/box[1][1];
                    G[2] = 2*M_PI*kz/box[2][2];
                    complex<double> SG_real = 0;
                    complex<double> SG_img = 0;
                    //real part
                    iter = omp_get_thread_num();
                    #pragma omp parallel for reduction(+:SG_real) schedule(dynamic) {
                        for(int i=0; i<n_atoms; i++){
                            double G_dot_r = G[0]*pos_ion[i][0] + G[1]*pos_ion[i][1] + G[2]*pos_ion[i][2];
                            SG_real += ion_charge[i]*cos(G_dot_r);
                        }
                    }
                    //img part
                    #pragma omp parallel for reduction(+:SG_img) schedule(dynamic){
                        for(int i=0; i<n_atoms; i++){
                            double G_dot_r = G[0]*pos_ion[i][0] + G[1]*pos_ion[i][1] + G[2]*pos_ion[i][2];
                            SG_img += ion_charge[i]*t*sin(G_dot_r);
                        }
                    }
                    SG = SG_real + SG_img;
                    double norm_SG = norm(SG);
                    double mod_G = sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]);
                    double damping_factor_exp = -(mod_G*mod_G)/(4*beta*beta);
                    reci_energy_h+= (exp(damping_factor_exp)*norm_SG)/(mod_G*mod_G);

                }
            }
        }
        //normalisation of energy
        reci_energy_h= reci_energy*2*M_PI/(box[0][0]*box[1][1]*box[2][2]);
        return reci_energy;
    }
#elif defined REDUCTION_KVECTOR
    