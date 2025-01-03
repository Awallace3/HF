#include "helper.hpp"
#include "input.hpp"
#include "omp.h"
#include "stdio.h"
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
// #include <Eigen/Dense>
#include <eigen3/Eigen/Dense>

using namespace std;

void HF_og() {

  // Specify Data Path
  /* std::string dataPath = "data/t1"; */
  std::string dataPath = "data/t0";
  double t1 = 1e-8, t2 = 1e-8;


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

  // Set Number of Electrons for a Neutral Molecule
  helper::getNumberOfElectrons(num_atoms, elements, &num_electrons);
  cout << "Number of Electrons: " << num_electrons << endl;

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

void HF(int num_atoms, double E = 0, double e_nuc = 0,
        std::vector<int> *elements = nullptr,
        std::vector<double> *eri = nullptr, Eigen::MatrixXd *coords = nullptr,
        Eigen::MatrixXd *T = nullptr, Eigen::MatrixXd *V = nullptr,
        Eigen::MatrixXd *S = nullptr

) {
  // Specify Data Path
  double t1 = 1e-8, t2 = 1e-8;

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

  cout << "Number of Electrons: " << num_electrons << endl;

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
  free(S);
  free(eri);
  printf("\nFreed Memory\n");
}

void timings(std::string dataPath, int num_threads) {
  int num_atoms;
  double E = 0, e_nuc = 0;
  time_t start, end;
  double itime, ftime, exec_time;

  std::vector<int> *elements = nullptr;
  std::vector<double> *eri = nullptr;

  Eigen::MatrixXd *coords = nullptr;
  Eigen::MatrixXd *T = nullptr;
  Eigen::MatrixXd *V = nullptr;
  Eigen::MatrixXd *S = nullptr;
  input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &S,
                    &e_nuc);
  cout << "S: " << endl << *S << endl;
  omp_set_num_threads(num_threads);
  Eigen::setNbThreads(num_threads);
  double totTime;
  omp_set_num_threads(num_threads);
  Eigen::setNbThreads(num_threads);
  start = clock();
  itime = omp_get_wtime();
  HF(num_atoms, E, e_nuc, elements, eri, coords, T, V, S);
  ftime = omp_get_wtime();
  end = clock();
  totTime = (double)(end - start);
  exec_time = ftime - itime;

  cout << "Time (CPU) : " << (double)(totTime / CLOCKS_PER_SEC) << endl;
  cout << "Time (USR) : " << exec_time << endl << endl ;
}

void timings_parrallel(std::string dataPath, int num_threads) {
  /* std::string dataPath = "data/t1"; */
  /* std::string dataPath = "data/t1"; */
  int num_atoms;
  double E = 0, e_nuc = 0;
  time_t start, end;
  double serial_t;
  double itime, ftime, exec_time;

  // Required code for which execution time needs to be computed

  std::vector<int> *elements = nullptr;
  std::vector<double> *eri = nullptr;

  Eigen::MatrixXd *coords = nullptr;
  Eigen::MatrixXd *T = nullptr;
  Eigen::MatrixXd *V = nullptr;
  Eigen::MatrixXd *S = nullptr;
  input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &S,
                    &e_nuc);
  // cout S
  omp_set_num_threads(1);
  Eigen::setNbThreads(1);
  start = clock();
  HF(num_atoms, E, e_nuc, elements, eri, coords, T, V, S);
  end = clock();
  serial_t = (double)(end - start) / CLOCKS_PER_SEC;
  cout << "Serial Time: " << serial_t << endl;

  input::gatherData(dataPath, num_atoms, &elements, &eri, &coords, &T, &V, &S,
                    &e_nuc);
  double omp_t;
  /* int num_threads = 10; */
  omp_set_num_threads(num_threads);
  Eigen::setNbThreads(num_threads);
  start = clock();
  itime = omp_get_wtime();
  HF(num_atoms, E, e_nuc, elements, eri, coords, T, V, S);
  ftime = omp_get_wtime();
  end = clock();
  omp_t = (double)(end - start);
  exec_time = ftime - itime;
  double ompSpeedUp = serial_t / exec_time;
  double eff = serial_t / (exec_time * num_threads);

  cout << "Serial Time (CPU) : " << serial_t << endl;
  cout << "OMP    Time (CPU) : " << (double)(omp_t / CLOCKS_PER_SEC) << endl;
  cout << "OMP    Time (USR) : " << exec_time << endl;
  cout << "Omp Speedup       : " << ompSpeedUp << endl;
  cout << "Parallel Efficieny: " << eff << endl;
}

int main(int argc, char *argv[]) {
  printf("\nRunning: %s\n\n", argv[0]);
  if (argc == 1) {
    printf("You must pass a data path and number of threads like "
           "below:\n\t./hf data/t1 4\n");
    return 1;
  }
  std::string dataPath = "";
  int numThreads = 1;
  if (argc >= 2) {
    dataPath = argv[1];
    numThreads = atoi(argv[2]);
  }
  cout << "Data Path: " << dataPath << endl;
  cout << "Num Threads: " << numThreads << endl;
  cout.precision(15);
  /* timings_parrallel(dataPath, numThreads); */
  // HF_og();
  timings(dataPath, numThreads);
  return 0;
}
