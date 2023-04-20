#ifndef HELPER_HPP
#define HELPER_HPP 1

#include "stdio.h"
#include <Eigen/Dense>
#include <vector>

namespace helper {
void orthoBasisSet(Eigen::MatrixXd *S, Eigen::MatrixXd *X);

void initialFockMatrix(Eigen::MatrixXd *S, Eigen::MatrixXd *H,
                       Eigen::MatrixXd *F);

void CMatrix(Eigen::MatrixXd *F, Eigen::MatrixXd *C);

void initialDensityMatrix(Eigen::MatrixXd *C, Eigen::MatrixXd *D, int num_electrons);

void initialEnergy(Eigen::MatrixXd *P, Eigen::MatrixXd *H, Eigen::MatrixXd *F,
                   double *E);

void getNumberOfElectrons(int num_atoms, std::vector<int> *elements,
                          int *num_electrons);

int indexIJKL(int i, int j, int k, int l);

void eriReducedCalc(std::vector<double> *eri, std::vector<double> *eriReduced);


} // namespace helper
#endif
