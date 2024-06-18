#include "libinclude.h"
#include "fundec.h"
// #include "omp.h"
#include "const.h"

#define PADDING_SIZE 8
// Uncomment only one of the following options
//#define USE_NAIVE 1
//#define USE_FALSE_SHARING_AND_PADDING 2
//#define USE_REDUCTION 3
#define USE_SYNCHRONIZATION_CONSTRUCT 4

#if defined USE_FALSE_SHARING_AND_PADDING
//* Avoiding false sharing through padding
double compute_real_energy(double **atomPositions, float *atomCharges, int numAtoms, double beta, float **simulationBox) {
    int totalThreads;
    double partialSums[NUM_THREADS][PADDING_SIZE] = {0};
    double totalEnergy = 0.0;
    omp_set_num_threads(thread::hardware_concurrency()); //Sets the number of OpenMP threads 
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num(); //Retrieves the ID of the current thread
        if (threadID == 0){ //master thread
            totalThreads = activeThreads;
        }

        for (int i = threadID, partialSum = 0; i < numAtoms; i += activeThreads) {
            for (int j = 0; j < i; j++) {
                if (i != j) {
                    double distance = dist(atomPositions, i, j, simulationBox);
                    partialSums[threadID][0] += (atomCharges[i] * atomCharges[j] * erfc(beta * distance)) / distance;
                }
            }
        }
    }

    for (int k = 0; k < totalThreads; k++) {
        totalEnergy += partialSums[k][0];
    }
    return totalEnergy;
}

#elif defined USE_SYNCHRONIZATION_CONSTRUCT
//* Synchronization construct with critical section
double compute_real_energy(double **atomPositions, float *atomCharges, int numAtoms, double beta, float **simulationBox) {
    int totalThreads;
    double totalEnergy = 0.0;
    omp_set_num_threads(thread::hardware_concurrency());

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        int activeThreads = omp_get_num_threads();
        if (threadID == 0) totalThreads = activeThreads;

        double localSum = 0.0;
        for (int i = threadID; i < numAtoms; i += activeThreads) {
            for (int j = 0; j < i; j++) {
                if (i != j) {
                    double distance = dist(atomPositions, i, j, simulationBox);
                    localSum += (atomCharges[i] * atomCharges[j] * erfc(beta * distance)) / distance;
                }
            }
        }

        #pragma omp critical
        totalEnergy += localSum;
    }
    return totalEnergy;
}

#elif defined USE_NAIVE
//* Single-threaded computation (naive approach)
double compute_real_energy(double **atomPositions, float *atomCharges, int numAtoms, double beta, float **simulationBox) {
    double totalEnergy = 0.0;

    for (int i = 0; i < numAtoms; i++) {
        for (int j = 0; j < i; j++) {
            if (i != j) {
                double distance = dist(atomPositions, i, j, simulationBox);
                totalEnergy += (atomCharges[i] * atomCharges[j] * erfc(beta * distance)) / distance;
            }
        }
    }
    return totalEnergy;
}

#elif defined USE_REDUCTION
//* Parallel reduction construct
double compute_real_energy(double **atomPositions, float *atomCharges, int numAtoms, double beta, float **simulationBox) {
    double totalEnergy = 0.0;
    omp_set_num_threads(thread::hardware_concurrency());

    #pragma omp parallel for simd schedule(runtime) reduction(+:totalEnergy)
    for (int i = 0; i < numAtoms; i++) {
        for (int j = 0; j < i; j++) {
            if (i != j) {
                double distance = dist(atomPositions, i, j, simulationBox);
                totalEnergy += (atomCharges[i] * atomCharges[j] * erfc(beta * distance)) / distance;
            }
        }
    }
    return totalEnergy;
}


#else
#error "Please define a parallelization method"
#endif