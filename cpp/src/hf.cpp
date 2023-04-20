#include "input.hpp"
#include "stdio.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {
  printf("Start\n");

  std::string filename = "data/water.xyz";
  /* std::vector<std::vector<double>> *coords; */

  std::vector<int> *elements = nullptr;
  std::vector<std::vector<double>> *coords = nullptr;
  int num_atoms;
  printf("Reading geometry...\n");
  input::readGeometry(filename, num_atoms, &elements, &coords);
  std::cout << "Number of atoms: " << num_atoms << std::endl;
  input::printElements(elements);
  input::print2dVector(coords);
  free(elements);
  printf("\nEnd\n");

  return 0;
}
