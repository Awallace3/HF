#ifndef INPUT_HPP
#define INPUT_HPP 1

#include "stdio.h"
#include <eigen3/Eigen/Dense>
#include <string>
#include <vector>
#include "helper.hpp"

namespace input {
void gatherData(std::string, int &, std::vector<int> **, std::vector<double> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **

);

void gatherData(std::string, int &, std::vector<int> **, std::vector<double> **,
                Eigen::MatrixXd **, Eigen::MatrixXd **, Eigen::MatrixXd **,
                Eigen::MatrixXd **, Eigen::MatrixXd **, double *

);

void gatherData(std::string, int &, std::vector<int> **, std::vector<double> **,
                Eigen::MatrixXd **, Eigen::MatrixXd **,
                Eigen::MatrixXd **, Eigen::MatrixXd **, double *

);

void readVector(std::string, std::vector<std::vector<double>> **);
void readVector(std::string, std::vector<double> **);
void readVector(std::string, Eigen::MatrixXd **);

/* void readERI(std::string, std::vector<double> **); */
void readERI(std::string fn, std::vector<double> **arr, int);

void readGeometry(std::string, int &, std::vector<int> **,
                  std::vector<std::vector<double>> **);
void readGeometry(std::string, int &, std::vector<int> **, Eigen::MatrixXd **);

void numAtoms(std::string, int &);

void printVector(std::vector<std::vector<double>> *);
void printVector(std::vector<double> *);
void printElements(std::vector<int> *);

void readNumber(std::string filename, double &number);
} // namespace input

#endif
