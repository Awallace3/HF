#include "input.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

void input::showHello() { std::cout << "Hello World!" << std::endl; }

void input::numAtoms(std::string filename, int &num_atoms) {
  std::ifstream file(filename);
  if (!file) {
    std::cout << "Could not open file " << filename << std::endl;
    return;
  }
  file >> num_atoms;
}

void input::readGeometry(std::string filename,
                         int &num_atoms,
                         std::vector<int> **elements,
                         std::vector<double> *x_coords,
                         std::vector<double> *y_coords,
                         std::vector<double> *z_coords
                         ) {
  std::ifstream file(filename);
  if (!file) {
    std::cout << "Could not open file " << filename << std::endl;
    return;
  }

  file >> num_atoms;
  *elements = new std::vector<int>(num_atoms);
  /* *x_coords = new std::vector<double>(num_atoms); */
  /* *y_coords = new std::vector<double>(num_atoms); */
  /* *z_coords = new std::vector<double>(num_atoms); */

  std::string line;
  std::getline(file, line); // read in the comment line

  for (int i = 0; i < num_atoms; ++i) {
    int el = -1;
    double x, y, z;
    file >> el >> x >> y >> z;

    /* elements.push_back(name); */
    (*elements)->at(i) = el;
    /* std::cout << elements->size() << std::endl; */
    /* x_coords.push_back(x); */
    /* y_coords.push_back(y); */
    /* z_coords.push_back(z); */
  }
  file.close();
  return;
}

void input::print2dVector(std::vector<std::vector<double>> *matrix) {
  for (int i = 0; u_int64_t(i) < matrix->size(); ++i) {
    for (int j = 0; u_int64_t(j) < matrix->at(i).size(); ++j) {
      std::cout << matrix->at(i).at(j) << " ";
    }
    std::cout << std::endl;
  }
}

void input::printElements(std::vector<int> *matrix) {
  for (int i = 0; u_int64_t(i) < matrix->size(); ++i) {
      std::cout << matrix->at(i) << " ";
    std::cout << std::endl;
  }
}

