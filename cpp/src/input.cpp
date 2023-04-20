#include "input.hpp"
#include "helper.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

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
      /* std::cout << value << " "; */
      j++;
    }
    /* std::cout << std::endl; */
    i++;
  }
  file.close();
  return;
}

void input::readVector(std::string fn, Eigen::MatrixXd **arr) {
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

  *arr = new Eigen::MatrixXd(count, count);

  int i = 0, j;
  while (getline(file, line)) {
    j = 0;
    std::stringstream ss(line);
    std::vector<double> values;
    double value;
    while (ss >> value) {
      /* std::cout << value << " " << i << " " << j << std::endl; */
      (**arr)(i, j) = value;
      j++;
    }
    /* std::cout << std::endl; */
    i++;
  }
  file.close();
  return;
}
void input::readVector(std::string fn, std::vector<double> **arr, int *num_electrons) {
    // TODO: need to read ERI into Nx4 matrix and use IJKL
    // indexing to build eri reduced matrix
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

void input::readERI(std::string fn, std::vector<double> **arr) {
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
  input::readERI(eriFN, eri);
  /* printVector(*eri); */
}

void input::gatherData(std::string dataPath, int &num_atoms,
                       std::vector<int> **elements,
                       std::vector<double> **eri,
                       Eigen::MatrixXd **coords,
                       Eigen::MatrixXd **T,
                       Eigen::MatrixXd **V,
                       Eigen::MatrixXd **e1,
                       Eigen::MatrixXd **overlap

) {
  std::string geom = dataPath + "/geom.xyz";
  std::string eriFN = dataPath + "/eri.dat";
  std::string TFN = dataPath + "/T.dat";
  std::string VFN = dataPath + "/V.dat";
  std::string e1FN = dataPath + "/e1.dat";
  std::string overlapFN = dataPath + "/overlap.dat";
  // Gathering Geometry
  input::readGeometry(geom, num_atoms, elements, coords);
  /* std::cout << "Number of atoms: " << num_atoms << std::endl; */
  /* input::printElements(*elements); */
  /* input::printVector(*coords); */

  // Gathering T, V, e1, overlap
  input::readVector(TFN, T);
  input::readVector(VFN, V);
  input::readVector(e1FN, e1);
  input::readVector(overlapFN, overlap);

  // Gathering eri
  /* input::readVector(eriFN, eri); */
  input::readERI(eriFN, eri);
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
                         Eigen::MatrixXd **coords) {
  std::ifstream file(filename);
  if (!file) {
    std::cout << "Could not open file " << filename << std::endl;
    return;
  }

  file >> num_atoms;
  *elements = new std::vector<int>(num_atoms);
  *coords = new Eigen::MatrixXd(num_atoms, 3);

  std::string line;
  std::getline(file, line); // read in the comment line

  for (int i = 0; i < num_atoms; ++i) {
    int el = -1;
    double x, y, z;
    file >> el >> x >> y >> z;
    (*elements)->at(i) = el;
    (**coords)(i, 0) = x;
    (**coords)(i, 1) = y;
    (**coords)(i, 2) = z;
  }
  file.close();
  return;
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
