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
  /* std::vector<std::vector<double>> *coords = nullptr; */
  /* std::vector<std::vector<double>> *T = nullptr; */
  /* std::vector<std::vector<double>> *V = nullptr; */
  /* std::vector<std::vector<double>> *e1 = nullptr; */
  /* std::vector<std::vector<double>> *S = nullptr; */

  Eigen::MatrixXd *coords = nullptr;
  Eigen::MatrixXd *T = nullptr;
  Eigen::MatrixXd *V = nullptr;
  Eigen::MatrixXd *e1 = nullptr;
  Eigen::MatrixXd *S = nullptr;
  // Read Data
  /* input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &e1, */
  /*                   &S); */
  input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &e1,
                    &S);
  printf("Read Data");
  // Starting HF Code
  // diagonalize S

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
