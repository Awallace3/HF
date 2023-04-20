#include "helper.hpp"
#include <Eigen/Dense>
#include <vector>

void helper::orthoBasisSet(Eigen::MatrixXd *S, Eigen::MatrixXd *X) {
  // Diagonalize S
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*S);
  Eigen::MatrixXd D = es.eigenvalues().asDiagonal();
  Eigen::MatrixXd U = es.eigenvectors();
  // Invert D
  for (int i = 0; i < D.rows(); i++) {
    D(i, i) = 1 / sqrt(D(i, i));
  }
  // Calculate X
  *X = U * D * U.transpose();
}

void helper::initialFockMatrix(Eigen::MatrixXd *S, Eigen::MatrixXd *H,
                       Eigen::MatrixXd *F) {
  // Calculate X
  Eigen::MatrixXd X;
  helper::orthoBasisSet(S, &X);
  // Calculate F
  *F = X.transpose() * *H * X;
}

void helper::CMatrix(Eigen::MatrixXd *F, Eigen::MatrixXd *C) {
  // Diagonalize F
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*F);
  Eigen::MatrixXd D = es.eigenvalues().asDiagonal();
  Eigen::MatrixXd U = es.eigenvectors();
  // Calculate C
  *C = U * D * U.transpose();
}

void helper::initialEnergy(Eigen::MatrixXd *P, Eigen::MatrixXd *H, Eigen::MatrixXd *F,
                   double *E) {
  // Calculate E
  *E = (0.5 * P->transpose() * (*H + *F)).sum();
}

void helper::getNumberOfElectrons(int num_atoms, std::vector<int> *elements,
                          int *num_electrons) {
  // Calculate number of electrons
  *num_electrons = 0;
  for (int i = 0; i < num_atoms; i++) {
    *num_electrons += elements->at(i);
  }
}

int helper::indexIJKL(int i, int j, int k, int l) {
  if (i > j)
    std::swap(i, j);
  if (k > l)
    std::swap(k, l);
  int ij = i * (i + 1) / 2 + j;
  int kl = k * (k + 1) / 2 + l;
  int ijkl = ij * (ij + 1) / 2 + kl;
  return ijkl;
}
