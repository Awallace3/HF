#include "input.hpp"
#include "stdio.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {
  printf("Start\n");

  std::string dataPath = "data";
  std::string geom = "data/geom.xyz";
  std::string eriFN = "data/eri.dat";
  std::string TFN = "data/T.dat";
  std::string VFN = "data/V.dat";
  std::string e1FN = "data/e1.dat";
  std::string overlapFN = "data/overlap.dat";

  int num_atoms;
  std::vector<int> *elements = nullptr;
  std::vector<double> *eri = nullptr;
  std::vector<std::vector<double>> *coords = nullptr;
  std::vector<std::vector<double>> *T = nullptr;
  std::vector<std::vector<double>> *V = nullptr;
  std::vector<std::vector<double>> *e1 = nullptr;
  std::vector<std::vector<double>> *S = nullptr;
  /* input::readGeometry(geom, num_atoms, &elements, &coords); */
  input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &e1,
                    &S);

  free(elements);
  free(coords);
  printf("\nEnd\n");
  return 0;
}
