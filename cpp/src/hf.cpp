#include "input.hpp"
#include "stdio.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "helper.hpp"

using namespace std;

int main() {
  // Specify Data Path
  std::string dataPath = "data";

  // Make pointers to store input data
  int num_atoms;
  std::vector<int> *elements = nullptr;
  std::vector<double> *eri = nullptr;

  Eigen::MatrixXd *coords = nullptr;
  Eigen::MatrixXd *T = nullptr;
  Eigen::MatrixXd *V = nullptr;
  Eigen::MatrixXd *e1 = nullptr;
  Eigen::MatrixXd *S = nullptr;

  // Read Data
  input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &e1,
                    &S);

  // Starting HF Code

  // Allocate Memory for Matrices
  Eigen::MatrixXd *H = nullptr;
  Eigen::MatrixXd *X = nullptr;
  Eigen::MatrixXd *F = nullptr;
  Eigen::MatrixXd *C = nullptr;
  Eigen::MatrixXd *D = nullptr;

  // Allocate memory for energy and electron count
  int num_electrons;
  double energy = 0;
  int eriSize;

  // Set Number of Electrons for a Neutral Molecule
  helper::getNumberOfElectrons(num_atoms, elements, &num_electrons);

  H = new Eigen::MatrixXd(T->rows(), T->cols());
  *H = *T + *V + *e1;
  /* cout << "H Matrix: " << endl << *H << endl; */

  // Orthogonalize Basis Set
  X = new Eigen::MatrixXd(S->rows(), S->cols());
  helper::orthoBasisSet(S, X);
  /* cout << "X Matrix: " << endl << *X << endl; */

  // Build initial Fock Matrix
  F = new Eigen::MatrixXd(H->rows(), H->cols());
  /* helper::initialFockMatrix(S, H, F); */
  *F = (*X).transpose() * *H * (*X);
  cout << "F Matrix: " << endl << *F << endl;

  C = new Eigen::MatrixXd(H->rows(), H->cols());
  helper::CMatrix(F, C);
  /* cout << "C Matrix: " << endl << *C << endl; */


  D = new Eigen::MatrixXd(H->rows(), H->cols());
  helper::initialDensityMatrix(C, D, num_electrons);
  cout << "D Matrix: " << endl << *D << endl;

  helper::initialEnergy(D, H, F, &energy);
  cout << "Energy: " << energy << endl;

  std::vector<double> *eri2 = nullptr;
  eriSize = H->rows() * H->cols() * H->rows() * H->cols();
  eri2 = new std::vector<double>(eriSize);
  /* helper::eriReducedCalc(eri, eri2); */





  // Free Allocations
  free(elements);
  free(coords);
  free(T);
  free(V);
  free(e1);
  free(S);
  free(eri);
  printf("\nFreed Memory\n");
  return 0;
}
