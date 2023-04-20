#include "input.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void input::readVector(std::string fn, std::vector<std::vector<double>> **arr) {
  std::ifstream file(fn);
  if (!file) {
    std::cout << "Could not open file: " << fn << std::endl;
    return;
  }
  std::string line;
  int count = 0;

  while (getline(file, line)) {
    count++;
  }
  file.clear();
  file.seekg(0);

  *arr = new std::vector<std::vector<double>>(count);
  for (int i = 0; i < count; ++i) {
    (*arr)->at(i) = std::vector<double>(count);
  }

  int i = 0, j;
  while (getline(file, line)) {
    j = 0;
    std::stringstream ss(line);
    std::vector<double> values;
    double value;
    while (ss >> value) {
      (*arr)->at(i).at(j) = value;
      std::cout << value << " ";
      j++;
    }
    std::cout << std::endl;
    i++;
  }
  file.close();
  return;
}

void input::readVector(std::string fn, std::vector<double> **arr) {
  std::ifstream file(fn);
  if (!file) {
    std::cout << "Could not open file: " << fn << std::endl;
    return;
  }
  std::string line;
  int count = 0;

  while (getline(file, line)) {
    count++;
  }
  file.clear();
  file.seekg(0);

  *arr = new std::vector<double>(count);

  int i = 0;
  while (getline(file, line)) {
    std::stringstream ss(line);
    ss >> (*arr)->at(i);
    i++;
  }
  file.close();
  return;
}

void input::gatherData(std::string dataPath, int &num_atoms,
                       std::vector<int> **elements, std::vector<double> **eri,
                       std::vector<std::vector<double>> **coords,
                       std::vector<std::vector<double>> **T,
                       std::vector<std::vector<double>> **V,
                       std::vector<std::vector<double>> **e1,
                       std::vector<std::vector<double>> **overlap

) {
  std::string geom = dataPath + "/geom.xyz";
  std::string eriFN = dataPath + "/eri.dat";
  std::string TFN = dataPath + "/T.dat";
  std::string VFN = dataPath + "/V.dat";
  std::string e1FN = dataPath + "/e1.dat";
  std::string overlapFN = dataPath + "/overlap.dat";
  // Gathering Geometry
  input::readGeometry(geom, num_atoms, elements, coords);
  std::cout << "Number of atoms: " << num_atoms << std::endl;
  input::printElements(*elements);
  input::printVector(*coords);

  // Gathering T, V, e1, overlap
  input::readVector(TFN, T);
  /* input::printVector(*T); */
  input::readVector(VFN, V);
  /* input::printVector(*V); */
  input::readVector(e1FN, e1);
  /* input::printVector(*e1); */
  input::readVector(overlapFN, overlap);
  input::readVector(eriFN, eri);
  printVector(*eri);
}

void input::numAtoms(std::string filename, int &num_atoms) {
  std::ifstream file(filename);
  if (!file) {
    std::cout << "Could not open file " << filename << std::endl;
    return;
  }
  file >> num_atoms;
}

void input::readGeometry(std::string filename, int &num_atoms,
                         std::vector<int> **elements,
                         std::vector<std::vector<double>> **coords) {
  std::ifstream file(filename);
  if (!file) {
    std::cout << "Could not open file " << filename << std::endl;
    return;
  }

  file >> num_atoms;
  *elements = new std::vector<int>(num_atoms);
  *coords = new std::vector<std::vector<double>>(num_atoms);
  for (int i = 0; i < num_atoms; ++i) {
    (*coords)->at(i) = std::vector<double>(3);
  }

  std::string line;
  std::getline(file, line); // read in the comment line

  for (int i = 0; i < num_atoms; ++i) {
    int el = -1;
    double x, y, z;
    file >> el >> x >> y >> z;
    (*elements)->at(i) = el;
    (*coords)->at(i).at(0) = x;
    (*coords)->at(i).at(1) = y;
    (*coords)->at(i).at(2) = z;
  }
  file.close();
  return;
}

void input::printVector(std::vector<std::vector<double>> *matrix) {
  std::cout << std::endl;
  for (int i = 0; u_int64_t(i) < matrix->size(); ++i) {
    for (int j = 0; u_int64_t(j) < matrix->at(i).size(); ++j) {
      std::cout.precision(12);
      std::cout << matrix->at(i).at(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void input::printVector(std::vector<double> *matrix) {
  std::cout << std::endl;
  for (int i = 0; u_int64_t(i) < matrix->size(); ++i) {
      std::cout.precision(12);
      std::cout << matrix->at(i) << " ";
  }
  std::cout << std::endl;
}

void input::printElements(std::vector<int> *matrix) {
  for (int i = 0; u_int64_t(i) < matrix->size(); ++i) {
    std::cout << matrix->at(i) << " ";
    std::cout << std::endl;
  }
}
