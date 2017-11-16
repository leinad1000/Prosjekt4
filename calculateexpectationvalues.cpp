#include "calculateexpectationvalues.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <armadillo>
#include "mpi.h"

using namespace std;
using namespace arma;
ofstream ofile;

// Initialize the seed and call the Mersienne algorithm
random_device rd;
mt19937_64 gen(rd());
uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

calculateExpectationValues::calculateExpectationValues(double E,
                                                       double M,
                                                       double* w,
                                                       double temperature,
                                                       int L,
                                                       int** spinMatrix,
                                                       long int mcs,
                                                       double* average,
                                                       int numberOfAcceptedConfigurations,
                                                       char* outfilename,
                                                       double initialTemperature,
                                                       double finalTemperature,
                                                       int mcs_start,
                                                       int RankProcess,
                                                       int NProcesses,
                                                       double* TotalExpectationValues
                                                       )
{

}

void calculateExpectationValues::initialE(int** spinMatrix,
                                          double& E,
                                          int L)
{
    E = 0.0;
    for (int x=0;x<L;x++) {
        for (int y=0;y<L;y++) {
            E -= (double) spinMatrix[x][y]*(spinMatrix[periodic(x,L,-1)][y] + spinMatrix[x][periodic(y,L,-1)]);
        }
    }
}

void calculateExpectationValues::initialM(int** spinMatrix,
                                          double& M,
                                          int L)
{
    M = 0.0;
    for (int x=0;x<L;x++) {
        for (int y=0;y<L;y++) {
        M += (double) spinMatrix[x][y];
        }
    }
}

int calculateExpectationValues::periodic(int spinNumber,
                                         int L,
                                         int add)
{
    return (int) ((spinNumber + L + add) % (L));
}

void calculateExpectationValues::metropolis(int L,
                                            int** spinMatrix,
                                            double& E,
                                            double& M,
                                            double* w,
                                            int& numberOfAcceptedConfigurations)
{
    for (int x=0;x<L;x++) {
        for (int y=0;y<L;y++) {
            // Finding random position in matrix
            int xPosition = floor(RandomNumberGenerator(gen)*L);
            int yPosition = floor(RandomNumberGenerator(gen)*L);
            // Finding deltaE
            int under = spinMatrix[xPosition][periodic(yPosition,L,-1)];
            int over = spinMatrix[xPosition][periodic(yPosition,L,1)];
            int venstre = spinMatrix[periodic(xPosition,L,-1)][yPosition];
            int hoyre = spinMatrix[periodic(xPosition,L,1)][yPosition];
            int deltaE = 2*spinMatrix[xPosition][yPosition]*(under + over + venstre + hoyre);

            // Updating M, E and spin at current position if test is positive

            if (RandomNumberGenerator(gen) <=w[deltaE + 8]) {
                spinMatrix[xPosition][yPosition] *= -1;
                numberOfAcceptedConfigurations += 1;
                M += (double) 2*spinMatrix[xPosition][yPosition];
                E += deltaE;

            }
        }
    }
}

void calculateExpectationValues::output(int L,
                                        long int mcs,
                                        double temperature,
                                        double* average,
                                        int numberOfAcceptedConfigurations)
{
    double norm = 1/((double) mcs);
    // Expectation values
    double averageE = average[0]*norm;
    double averageEsquared = average[1]*norm;
    double averageM = average[2]*norm;
    double averageMsquared = average[3]*norm;
    double averageAbsoluteM = average[4]*norm;
    double C_v = (1.0/(temperature*temperature))*(averageEsquared - averageE*averageE);
    double chi = (1.0/temperature)*(averageMsquared - averageAbsoluteM*averageAbsoluteM);
    // Variances
    //double varianceE = (averageEsquared - pow(averageE,2.0))/(L*L);
    //double varianceM = (averageMsquared - pow(averageM,2.0))/(L*L);
    // Writing to file
    //L = 1;
    ofile << setw(15) << setprecision(8) << mcs << "\t"
          << setw(15) << setprecision(8) << numberOfAcceptedConfigurations << "\t"
          << setw(15) << setprecision(8) << averageE/(L*L) << "\t"
          << setw(15) << setprecision(8) << averageAbsoluteM/(L*L) << "\t"
          << setw(15) << setprecision(8) << temperature << "\t"
          << setw(15) << setprecision(8) << C_v/(L*L) << "\t"
          << setw(15) << setprecision(8) << chi/(L*L) << endl;
}

void calculateExpectationValues::doCalculation(double E,
                                               double M,
                                               double* w,
                                               double temperature,
                                               int L,
                                               int** spinMatrix,
                                               long int mcs,
                                               double* average,
                                               int numberOfAcceptedConfigurations,
                                               char* outfilename)
{

    // Opening object for writing to file
    ofile.open(outfilename);

    // Resetting for each temperature
    initialE(spinMatrix, E, L);
    initialM(spinMatrix, M, L);

    for (int i=0;i<17;i++) {
        w[i]=0;
    }
    for (int i=0;i<17;i+=4) {
        w[i] = exp(-(i-8)/temperature);
    }
    for (int i=0;i<5;i++) {
        average[i] = 0.0;
    }

    // Looping over Monte Carlo Cycles
    for (long int i=1;i<=mcs;i++) {
        metropolis(L, spinMatrix, E, M, w, numberOfAcceptedConfigurations); // Producing E and M
        average[0] += E;
        average[1] += E*E;
        average[2] += M;
        average[3] += M*M;
        average[4] += fabs(M);

        output(L, i, temperature, average, numberOfAcceptedConfigurations); // Printing to file for question c)
    }

    //output(L, mcs, temperature, average, numberOfAcceptedConfigurations); // Printing to file for question b)

    // Closing outputfile
    ofile.close();
}

void calculateExpectationValues::question_d(double E,
                                               double M,
                                               double* w,
                                               double temperature,
                                               int L,
                                               int** spinMatrix,
                                               long int mcs,
                                               double* average,
                                               int numberOfAcceptedConfigurations,
                                               char* outfilename)
{
    // Opening object for writing to file
    ofile.open(outfilename);

    // Resetting for each temperature
    initialE(spinMatrix, E, L);
    initialM(spinMatrix, M, L);

    for (int i=0;i<17;i++) {
        w[i]=0;
    }
    for (int i=0;i<17;i+=4) {
        w[i] = exp(-(i-8)/temperature);
    }
    for (int i=0;i<5;i++) {
        average[i] = 0.0;
    }

    int NumberOfOutcomes = 0;
    int N = 401; // Number of possible outcomes; from -800 to 800 in steps of 4
    int indeks;
    int *countingOutcomes = new int[N];
    int *possibleOutcomes = new int[N];
    double *P = new double[N]; // Probability
    for (int i=0;i<N;i++) {
        countingOutcomes[i] = 0;
        possibleOutcomes[i] = -800  + i*4;
    }

    // Looping over Monte Carlo Cycles
    for (long int i=1;i<=mcs;i++) {
        metropolis(L, spinMatrix, E, M, w, numberOfAcceptedConfigurations); // Producing E and M
        average[0] += E;
        average[1] += E*E;
        average[2] += M;
        average[3] += M*M;
        average[4] += fabs(M);

        if (i>10000) {
            indeks = (E + 800)/4;
            countingOutcomes[indeks] += 1;
            NumberOfOutcomes += 1;
        }

    }

    for (int i=0;i<N;i++) {
        P[i] = ((double) countingOutcomes[i])/((double) NumberOfOutcomes);
        ofile << setw(15) << setprecision(8) << possibleOutcomes[i] << "\t"
              << setw(15) << setprecision(8) << P[i] << endl;

    }

    delete [] countingOutcomes;
    delete [] possibleOutcomes;
    delete [] P;

    // Finding standard-deviation
    double norm = 1/((double) mcs);
    double averageE = average[0]*norm;
    double averageEsquared = average[1]*norm;
    double varianceE = (averageEsquared - averageE*averageE);
    cout << "Standardavviket til E =" << " "<<sqrt(varianceE) << endl;

    // Printet ut for mcs=1e6: Standardavviket til E = 56.9772

    // Closing outputfile
    ofile.close();
}

void calculateExpectationValues::question_e(double E,
                                            double M,
                                            double* w,
                                            double temperature,
                                            int L,
                                            int** spinMatrix,
                                            long int mcs,
                                            double* average,
                                            int numberOfAcceptedConfigurations,
                                            char* outfilename,
                                            double initialTemperature,
                                            double finalTemperature,
                                            int mcs_start,
                                            int RankProcess,
                                            int NProcesses,
                                            double* TotalExpectationValues)
{

    if (RankProcess == 0 && temperature == initialTemperature) {
        ofile.open(outfilename);
    }

    // Resetting for each temperature
    initialE(spinMatrix, E, L);
    initialM(spinMatrix, M, L);

    for (int i=0;i<17;i++) {
        w[i]=0;
    }
    for (int i=0;i<17;i+=4) {
        w[i] = exp(-(i-8)/temperature);
    }
    for (int i=0;i<5;i++) {
        average[i] = 0.0;
        TotalExpectationValues[i] = 0.0;
    }

    // Looping over Monte Carlo Cycles
    for (long int i=1;i<=mcs;i++) {
        metropolis(L, spinMatrix, E, M, w, numberOfAcceptedConfigurations); // Producing E and M
        if (i>mcs_start) {
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);

        }
    }

    for( int i =0; i < 5; i++){
      MPI_Reduce(&average[i], &TotalExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    for (int i=0;i<5;i++) {
        TotalExpectationValues[i] = TotalExpectationValues[i]/((double) NProcesses);
    }

    if ( RankProcess == 0) {
        output(L, mcs-mcs_start, temperature, TotalExpectationValues, numberOfAcceptedConfigurations);
    }

    if (RankProcess == 0 && temperature == finalTemperature) {
        ofile.close();
    }
}

