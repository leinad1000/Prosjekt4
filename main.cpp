#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include "random"
#include <stdlib.h>
#include <random>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <armadillo>
#include <cmath>
#include <string>
#include "mpi.h"
#include "calculateexpectationvalues.h"

using namespace std;
using namespace arma;

// Declaring function
void fillUpMatrix(int, bool, int**);

int main(int argc, char *argv[])
{

    // Variables
    int L;
    long int numberOfAcceptedConfigurations;
    long int mcs;
    double temperature, E, M;
    bool RandomStart = true;
    bool notRandomStart = false;
    char* outfilename;

    L = 20;
    temperature = 2.4;
    mcs = 1e6;

    E = M = 0.0;
    numberOfAcceptedConfigurations = 0;

    double *w = new double[17];
    double *average = new double[5];
    double *TotalExpectationValues = new double[5];

    // Reading in outputfile. Abort if there are too few command line arguments
    if (argc<2){
        cout << "Bad usage: " << argv[0] << "read also outputfile on the same line"<<endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }

    // Main part
    calculateExpectationValues CEV;

    
    // Use this for questions b) and c)
    //Declaring arrays and matrices
    int **spinMatrix = new int*[L];
    for(int i = 0; i < L; i++){
        spinMatrix[i] = new int[L];
    }

    fillUpMatrix(L, notRandomStart, spinMatrix);
    CEV.doCalculation(E, M, w, temperature, L, spinMatrix, mcs, average, numberOfAcceptedConfigurations, outfilename);

    /*
    // Test for question b)
    double Z = 12 + 2*exp(-8.0) + 2*exp(8.0);
    double excactE = (16*exp(-8.0) - 16*exp(8.0))/Z;
    double excactAbsoluteM = (8*exp(8) + 16)/Z;
    double excactC_v = (1/(temperature*temperature))*((384*(exp(-8/temperature)+exp(8/temperature))+256)/(pow(6+exp(-8/temperature)+exp(8/temperature),2.0)));
    double excactChi = (1/temperature)*((192*exp(8/temperature)+64*exp(-8/temperature) + 192)/pow((12+2*exp(-8/temperature)+2*exp(8/temperature)),2.0));
    cout << "E = "  << excactE << endl;
    cout << "M = "  << excactAbsoluteM << endl;
    cout << "C_v =" << excactC_v << endl;
    cout << "Chi =" << excactChi << endl;
    */

    // Use this for question d)
    //CEV.question_d(E, M, w, temperature, L, spinMatrix, mcs, average, numberOfAcceptedConfigurations, outfilename);

    /*
    // Question e)

    int NProcesses, RankProcess;
    int mcs_start;
    double initialTemperature, finalTemperature, epsilon, increment;

    //  MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &NProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &RankProcess);

   if (RankProcess == 0 && argc < 2) {
        cout << "Bad Usage: " << argv[0] <<
          " read in output file" << endl;
        exit(1);
    }
    if (RankProcess == 0 && argc == 2) {
        outfilename = argv[1];
        L = 60;
        mcs = 2*1e5;
        initialTemperature = 2.200;
        finalTemperature = 2.362;
        increment = 0.0030;
        epsilon = 0.00005;
        mcs_start = 30000;
    }

    // broadcast to all nodes common variables since only master node reads from command line
    MPI_Bcast (&mcs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initialTemperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&finalTemperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&increment, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&mcs_start, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Declaring arrays and matrices
    int **spinMatrix = new int*[L];
    for(int i = 0; i < L; i++){
        spinMatrix[i] = new int[L];
    }

    for (double temperature=initialTemperature;temperature<finalTemperature+epsilon;temperature+=increment) {
        fillUpMatrix(L, RandomStart, spinMatrix);
        CEV.question_e(E, M, w, temperature, L, spinMatrix, mcs, average, numberOfAcceptedConfigurations, outfilename,initialTemperature, finalTemperature, mcs_start, RankProcess, NProcesses, TotalExpectationValues);
    }

    // End MPI
    MPI_Finalize ();
    */

    // Deleting arrays and matrices
    delete [] w;
    delete [] average;
    delete [] TotalExpectationValues;

    for(int i = 0; i < L; ++i) {
        delete [] spinMatrix[i];
    }

    delete [] spinMatrix;

    return 0;
}


void fillUpMatrix(int L, bool random, int** spinMatrix) {
    // Initialize the seed and call the Mersienne algorithm
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    // All spins pointing up
    if (!random) {
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                spinMatrix[i][j] = 1;
            }
        }
    }

    // Random spins
    else {
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                if (RandomNumberGenerator(gen) < 0.5) {
                    spinMatrix[i][j] = -1;
                }
                else {
                    spinMatrix[i][j] = 1;
                }
            }
        }
    }
}

