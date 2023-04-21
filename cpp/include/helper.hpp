#ifndef HELPER_HPP
#define HELPER_HPP 1

#include "stdio.h"
#include <Eigen/Dense>
#include <vector>

namespace helper {
void orthoS(Eigen::MatrixXd *S, Eigen::MatrixXd *X);

void initialFockMatrix(Eigen::MatrixXd *S, Eigen::MatrixXd *H,
                       Eigen::MatrixXd *F);

void getC_0_prime(Eigen::MatrixXd *F, Eigen::MatrixXd *C);

void updateDensityMatrix(Eigen::MatrixXd *C, Eigen::MatrixXd *D,
                          int num_electrons);

void computeEnergy(Eigen::MatrixXd *D, Eigen::MatrixXd *H, Eigen::MatrixXd *F,
                   double *E);

void getNumberOfElectrons(int num_atoms, std::vector<int> *elements,
                          int *num_electrons);

int indexIJKL(int i, int j, int k, int l);

void eriReducedCalc(std::vector<double> *eri, std::vector<double> *eriReduced);

void updateFockMatrix(Eigen::MatrixXd *H,
        Eigen::MatrixXd *D,
        Eigen::MatrixXd *F,
        std::vector<double> *eri);

void SCF(std::vector<double> *eri,
                 Eigen::MatrixXd *S_12, Eigen::MatrixXd *H, Eigen::MatrixXd *F,
                 Eigen::MatrixXd *C, Eigen::MatrixXd *D,
                 Eigen::MatrixXd *C_0_prime,
                 int num_electrons,
                 double *E, double e_nuc, double t1, double t2);

} // namespace helper
#endif
