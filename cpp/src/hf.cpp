#include "input.hpp"
#include "stdio.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>

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
  // diagonalize S
  S->diagonalSize();

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
