---
title: HF Final Project
date: 2023-04-24
author: Austin M. Wallace
reference-section-title: "References"
link-citations: true
link-bibliography: true
header-includes:
- \usepackage{cancel}
- \usepackage{fancyvrb}
- \usepackage{listings}
- \usepackage{amsmath}
- \usepackage[version=4]{mhchem}
- \usepackage{xcolor}
- \lstset{breaklines=true}
- \lstset{language=[Motorola68k]Assembler}
- \lstset{basicstyle=\small\ttfamily}
- \lstset{extendedchars=true}
- \lstset{tabsize=2}
- \lstset{columns=fixed}
- \lstset{showstringspaces=false}
- \lstset{frame=trbl}
- \lstset{frameround=tttt}
- \lstset{framesep=4pt}
- \lstset{numbers=left}
- \lstset{numberstyle=\tiny\ttfamily}
- \lstset{postbreak=\raisebox{0ex}[0ex][0ex]{\ensuremath{\color{red}\hookrightarrow\space}}}
- \lstset{keywordstyle=\color[rgb]{0.13,0.29,0.53}\bfseries}
- \lstset{stringstyle=\color[rgb]{0.31,0.60,0.02}}
- \lstset{commentstyle=\color[rgb]{0.56,0.35,0.01}\itshape}
---

# Hartree-Fock Overview
The objective of Hartree-Fock (HF) is to solve the Schrodinger Equation through
the usage of optimizing coefficients for a linear combination of atomic
orbitals (LCAO) to create molecular orbitals that minimize the electronic
energy. Because this procedure depends on the electronic energy, the nuclear
energy is not added until the very end of the self-consistent field procedure
is complete. Within this procedure, orbitals can be optimized by minimizing the
electronic energy due to the varitaional theorem.

Effectively, this process relies on nuclear coordinates, charge, multiplicity,
and a basis set -- the present work uses atom-centered Gaussian functions to
describe atomic orbitals. From here, the overalp of atomic orbitals is computed
to forumlate S along with compiting kinetic (T) and potential (V) energy
intregals at the beginning. Next, the Core Hamiltonian is formed $H_{\mu\nu} =
T_{\mu\nu} + V_{\mu\nu}$. Before constructing the core Fock matrix (Equation 1) for a naive initial density matrix guess (Equation 3),
one must orthogonalize $S^{-1/2}$ through diagonalizing S and transforming 
to $S^{-1/2}$.
\begin{equation}
F_0' = (S^{-1/2})^\dagger H S^{-1/2}
\end{equation}
Through diagonlizing the core Fock matrix with the orthogonalized basis, one
can construct $C_0$.
\begin{equation}
C_0 = S^{-1/2}C_0'
\end{equation}

Using these sorted eigenvectors, we compute our density guess.
\begin{equation}
D_{\mu\nu} = \sum_{i}^{N/2} C_{\mu i} C_{\nu i}
\end{equation}
With these inital values, the SCF procedure may begin, where we iteratively compute a new
Fock matrix (Equation 4), compute electronic energy (Equation 5), transform
Fock matrix to the orthonormal basis (Equation 6), diagonlaize F, (Equation 7) update C (Equation 8), and
update the density (Equation 9) until converging the energy.
\begin{equation}
F_{\mu\nu} = H_{\mu\nu} + \sum_{\rho\sigma}^{AO}D_{\rho\sigma}{2 [\mu\nu| \rho\sigma] - [\mu\rho | \nu\sigma]}
\end{equation}
\begin{equation}
E = \sum_{\mu\nu}^{AO} D_{\mu\nu} (H_{\mu\nu} + F_{\mu\nu})
\end{equation}
\begin{equation}
F' = (S^{-1/2})^\dagger F S^{-1/2}
\end{equation}
\begin{equation}
C'^\dagger F' C' = \epsilon
\end{equation}
\begin{equation}
C = S^{-1/2} C'
\end{equation}
\begin{equation}
D_{\mu\nu} = \sum_{i}^{N/2} C_{\mu i} C_{\nu i}
\end{equation}
At the end of each iteration, convergence is checked by computing the energy
and comparing with the previous iteration and comparing with a user defined
threshold to break the iterative process. Once below the threshold, the
electronic energy is added to the nuclear energy to get the HF energy.

# Specific Implementation Details
I have implemented the abstract described HF procedure above in C++ through the
usage of Eigen3 for performing linear algebra operations with Lapack under the
hood. C++ was selected due to having strong support with OpenMP and MPI for
parallelization and to improve fundamental knowledege for aspring to become a
Psi4 developer.

The project allowed me to learn how to create my own cmake files for building
the library across OSX and Ubuntu with ease using the
[pitchfork](https://api.csswg.org/bikeshed/?force=1&url=https://raw.githubusercontent.com/vector-of-bool/pitchfork/develop/data/spec.bs)
convention. Using CMakeLists.txt took some initial learning but provides a
seamless path for finding libraries and repeatedly building the project during
development. Additionally, it makes it easier for others to know what is needed
for getting the project to run on their own machine. 

The initial test case for getting the library working follows the "Coding
Strategy #1" from canvas; however, the T, V, S, and ERI sections were broken
into separate files for a simplier parsing function to generate Eigen::MatrixXd
objects. The Eigen library provided an easy option for solving for eigenvalues
and eigenvectors through passing your matrix to the SelfAdjointEigenSolver
object initialization. The eigenvectors and eigenvalues already returned sorted
and have methods to reverse the sorting if needed.

# Code
The project is available on
[GitHub](https://github.com/Awallace3/HF/tree/master/cpp),
where I plan to continue to implement parallelization
towards a distributed HF code.

The most important files are `hf.cpp`, `input.cpp`, and
`helper.cpp` which are displayed below.

```cpp
#include "helper.hpp"
#include "input.hpp"
#include "omp.h"
#include "stdio.h"
#include <Eigen/Dense>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void serial() {
  // Specify Data Path
  /* std::string dataPath = "data/t1"; */
  std::string dataPath = "data/t0";
  double t1 = 1e-8, t2 = 1e-8;

  //

  // Make pointers to store input data
  int num_atoms;
  double E = 0, e_nuc = 0;
  std::vector<int> *elements = nullptr;
  std::vector<double> *eri = nullptr;

  Eigen::MatrixXd *coords = nullptr;
  Eigen::MatrixXd *T = nullptr;
  Eigen::MatrixXd *V = nullptr;
  Eigen::MatrixXd *e1 = nullptr;
  Eigen::MatrixXd *S = nullptr;

  // Read Data
  input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &e1,
                    &S, &e_nuc);
  cout << "e_nuc: " << e_nuc << endl;

  // Starting HF Code

  // Allocate Memory for Matrices
  Eigen::MatrixXd *H = nullptr;
  Eigen::MatrixXd *S_12 = nullptr;
  Eigen::MatrixXd *F = nullptr;
  Eigen::MatrixXd *C = nullptr;
  Eigen::MatrixXd *C_0_prime = nullptr;
  Eigen::MatrixXd *D = nullptr;

  // Allocate memory for energy and electron count
  int num_electrons;
  /* int eriSize; */

  // Set Number of Electrons for a Neutral Molecule
  helper::getNumberOfElectrons(num_atoms, elements, &num_electrons);

  H = new Eigen::MatrixXd(T->rows(), T->cols());
  *H = *T + *V;
  cout << endl << "H Matrix: " << endl << endl << *H << endl;

  // Orthogonalization of S
  cout << endl << "S Matrix: " << endl << endl << *S << endl;
  S_12 = new Eigen::MatrixXd(S->rows(), S->cols());
  helper::orthoS(S, S_12);
  cout << endl << "S_12 Matrix: " << endl << endl << *S_12 << endl;

  // Build initial Fock Matrix
  F = new Eigen::MatrixXd(H->rows(), H->cols());
  *F = (*S_12).transpose() * *H * (*S_12);
  cout << endl << "F Matrix: " << endl << endl << *F << endl;

  C_0_prime = new Eigen::MatrixXd(H->rows(), H->cols());
  helper::getC_0_prime(F, C_0_prime);
  cout << endl << "C_0_prime Matrix: " << endl << endl << *C_0_prime << endl;

  C = new Eigen::MatrixXd(H->rows(), H->cols());
  *C = (*S_12) * (*C_0_prime);
  cout << endl << "C Matrix: " << endl << endl << *C << endl;

  D = new Eigen::MatrixXd(H->rows(), H->cols());
  helper::updateDensityMatrix(C, D, num_electrons);
  cout << endl << "D Matrix: " << endl << endl << *D << endl;

  helper::SCF(eri, S_12, H, F, C, D, C_0_prime, num_electrons, &E, e_nuc, t1,
              t2);
  E += e_nuc;
  cout.precision(15);
  cout << endl << "Final HF Energy: " << E << endl;

  // Final HF Energy:   -74.9659010585405
  // Total Energy =     -74.9659900701199433

  // Free Allocations
  free(elements);
  free(coords);
  free(T);
  free(V);
  free(e1);
  free(S);
  free(eri);
  printf("\nFreed Memory\n");
}

int main() {
  omp_set_num_threads(1);
  Eigen::setNbThreads(1);
  time_t start, end;
  double serial_t;
  start = clock();
  serial();
  end = clock();
  serial_t = (double)(end - start);
  cout << "Serial Time: " << (serial_t / CLOCKS_PER_SEC) << endl;
  double omp_t;
  int num_threads = 2;
  omp_set_num_threads(num_threads);
  Eigen::setNbThreads(num_threads);
  start = clock();
  serial();
  end = clock();
  omp_t = (double)(end - start);
  cout << "Omp Time: " << (double) (omp_t / CLOCKS_PER_SEC) << endl;
  cout << "Omp Speedup: " << (double)(serial_t / omp_t)
       << endl;

  return 0;
}
```

```cpp
#include "input.hpp"
#include "helper.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

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
void input::readVector(std::string fn, std::vector<double> **arr) {
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

void input::readERI(std::string fn, std::vector<double> **arr, int n_basis) {
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
  int arrSize =
      /* n_basis * (n_basis + 1) / 2  * (n_basis + 2) / 2 * (n_basis  + 3) / 2 ;
       */
      /* n_basis * n_basis * n_basis * n_basis / 8; */
      /* helper::indexIJKL(n_basis, n_basis, n_basis, n_basis); */
      helper::indexIJKL(n_basis - 1, n_basis - 1, n_basis - 1, n_basis - 1) + 1;

  cout << "arrSize: " << arrSize << endl;
  cout << "count: " << count << endl;

  *arr = new std::vector<double>(count);

  int i, j, k, l;
  double value;
  int ijkl;
  while (getline(file, line)) {
    std::stringstream ss(line);
    ss >> i >> j >> k >> l >> value;
    ijkl = helper::indexIJKL(i, j, k, l);
    /* std::cout << ijkl << " " << i << " " << j << " " << k << " " << l << " "
     */
    /*           << value << std::endl; */
    (*arr)->at(ijkl) = value;
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
  /* printVector(*eri); */
}

void input::gatherData(std::string dataPath, int &num_atoms,
                       std::vector<int> **elements, std::vector<double> **eri,
                       Eigen::MatrixXd **coords, Eigen::MatrixXd **T,
                       Eigen::MatrixXd **V, Eigen::MatrixXd **e1,
                       Eigen::MatrixXd **overlap, double *enuc) {
  std::string geom = dataPath + "/geom.xyz";
  std::string eriFN = dataPath + "/eri.dat";
  std::string TFN = dataPath + "/T.dat";
  std::string VFN = dataPath + "/V.dat";
  std::string e1FN = dataPath + "/e1.dat";
  std::string overlapFN = dataPath + "/overlap.dat";
  std::string enucFN = dataPath + "/enuc.dat";
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
  int n_basis = (*T)->rows();
  cout << "n_basis: " << n_basis << endl;
  input::readERI(eriFN, eri, n_basis);
  input::readNumber(enucFN, *enuc);
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

void input::readNumber(std::string filename, double &number) {
  std::ifstream file(filename);
  if (!file) {
    std::cout << "Could not open file " << filename << std::endl;
    return;
  }
  file >> number;
  file.close();
}
```

```cpp
#include "helper.hpp"
#include "stdio.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

void helper::orthoS(Eigen::MatrixXd *S, Eigen::MatrixXd *S12) {
  // Diagonalize S
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*S);
  /* Eigen::MatrixXd D = es.eigenvalues().asDiagonal(); */
  Eigen::MatrixXd LAMBDA = es.eigenvalues().asDiagonal();
  Eigen::MatrixXd U = es.eigenvectors();
  // Invert D
  for (int i = 0; i < LAMBDA.rows(); i++) {
    LAMBDA(i, i) = 1 / sqrt(LAMBDA(i, i));
  }
  // Calculate X
  *S12 = U * LAMBDA * U.transpose();
}

void helper::initialFockMatrix(Eigen::MatrixXd *X, Eigen::MatrixXd *H,
                               Eigen::MatrixXd *F) {
  // Calculate F
  *F = (*X).transpose() * *H * (*X);
}

void helper::getC_0_prime(Eigen::MatrixXd *F, Eigen::MatrixXd *C) {
  // Diagonalize F
  SelfAdjointEigenSolver<MatrixXd> eigensolver(*F);
  if (eigensolver.info() != Success)
    abort(); // check for errors
  *C = eigensolver.eigenvectors();
}

void helper::computeEnergy(Eigen::MatrixXd *D, Eigen::MatrixXd *H,
                           Eigen::MatrixXd *F, double *E) {
  // Calculate E
  // TODO: fix this
  *E = 0;
  /* for (int i = 0; i < H->rows(); i++) { */
  /*   for (int j = 0; j < H->rows(); j++) { */
  /*     *E += (*D)(i, j) * ((*H)(i, j) + (*F)(i, j)); */
  /*   } */
  /* } */
  *E = (*D).cwiseProduct((*H) + (*F)).sum();
  /* *E = (*D * ( (*H) + (*E))); */
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
  if (j > i){
    std::swap(i, j);
  }
  if (l > k){
    std::swap(k, l);
  }
  int ij = i * (i + 1) / 2 + j;
  int kl = k * (k + 1) / 2 + l;
  if (ij < kl){
    std::swap(ij, kl);
  }
  int ijkl = ij * (ij + 1) / 2 + kl;
  return ijkl;
}

void helper::updateDensityMatrix(Eigen::MatrixXd *C, Eigen::MatrixXd *D,
                                 int num_electrons) {
  // Calculate D
    for (int i = 0; i < C->rows(); i++) {
      for (int j = 0; j < C->rows(); j++) {
        (*D)(i, j) = 0;
        for (int k = 0; k < num_electrons / 2; k++) {
          (*D)(i, j) += (*C)(i, k) * (*C)(j, k);
        }
      }
    }
}

// TODO: finish this function
void helper::eriReducedCalc(std::vector<double> *eri,
                            std::vector<double> *eriReduced) {
  /* for (int i =0; i < eri->size(); i++){ */
  /*     eriReduced->at(eri->at(i)); */
  /* } */
}

void helper::updateFockMatrix(Eigen::MatrixXd *H, Eigen::MatrixXd *D,
                              Eigen::MatrixXd *F, std::vector<double> *eri) {
  // Update Fock Matrix
  *F = *H;
  for (int mu = 0; mu < H->rows(); mu++) {
    for (int nu = 0; nu < H->cols(); nu++) {
      for (int rho = 0; rho < H->rows(); rho++) {
        for (int sig = 0; sig < H->cols(); sig++) {
          (*F)(mu, nu) += (*D)(rho, sig) *
                          (2 * eri->at(helper::indexIJKL(mu, nu, rho, sig)) -
                           eri->at(helper::indexIJKL(mu, rho, nu, sig)));
        }
      }
    }
  }
}

void helper::SCF(std::vector<double> *eri, Eigen::MatrixXd *S_12,
                 Eigen::MatrixXd *H, Eigen::MatrixXd *F, Eigen::MatrixXd *C,
                 Eigen::MatrixXd *D, Eigen::MatrixXd *C_0_prime,
                 int num_electrons, double *E, double e_nuc, double t1,
                 double t2) {
  // Calculate SCF

  bool converged = false;
  double E2 = 0;
  int iter = 0, max_iter = 100;
  while (!converged) {
    // Update Fock Matrix
    helper::updateFockMatrix(H, D, F, eri);
    /* cout << endl <<"F Matrix: " << endl << endl <<*F << endl; */
    /* cout << endl <<"E: " << endl << endl <<*E << endl; */
    helper::computeEnergy(D, H, F, E);
    *F = (*S_12).transpose() * *F * (*S_12);

    helper::getC_0_prime(F, C_0_prime);
    /* cout << endl <<"C_0_prime Matrix: " << endl <<endl << *C_0_prime << endl;
     */

    *C = (*S_12) * (*C_0_prime);
    /* cout << endl <<"C Matrix: " << endl <<endl << *C << endl; */
    helper::updateDensityMatrix(C, D, num_electrons);
    /* cout << endl <<"D Matrix: " << endl << endl <<*D << endl; */

    cout << "iter: " << iter << " Energy: " << *E << " Delta E: " << (*E - E2) << endl;
    if (abs(*E - E2) < t1) {
      converged = true;
    } else if (iter > max_iter) {
      cout << "Max iterations reached" << endl;
      converged = true;
    } else {
      E2 = *E;
    }
    iter++;
    /* converged = true; */
  }
}
```








