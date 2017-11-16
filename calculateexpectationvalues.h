#ifndef CALCULATEEXPECTATIONVALUES_H
#define CALCULATEEXPECTATIONVALUES_H


class calculateExpectationValues
{
public:
    calculateExpectationValues(double E,
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
                               );

    calculateExpectationValues() {}

    void doCalculation (double E,
                        double M,
                        double* w,
                        double temperature,
                        int L,
                        int** spinMatrix,
                        long int mcs,
                        double* average,
                        int numberOfAcceptedConfigurations,
                        char* outfilename);
    void question_d    (double E,
                        double M,
                        double* w,
                        double temperature,
                        int L,
                        int** spinMatrix,
                        long int mcs,
                        double* average,
                        int numberOfAcceptedConfigurations,
                        char* outfilename);
    void question_e    (double E,
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
                        double* TotalExpectationValues);
    void output        (int L,
                        long int mcs,
                        double temperature,
                        double* average,
                        int numberOfAcceptedConfigurations);
    void metropolis    (int L,
                        int** spinMatrix,
                        double& E,
                        double& M,
                        double * w,
                        int& numberOfAcceptedConfigurations);
    void initialE      (int** spinMatrix,
                        double& E,
                        int L);
    void initialM      (int** spinMatrix,
                        double& M,
                        int L);
    int periodic       (int spinNumber,
                        int L,
                        int add);
};

#endif // CALCULATEEXPECTATIONVALUES_H
