#ifndef INPUT_HPP
#define INPUT_HPP 1

#include "stdio.h"
#include <string>
#include <vector>
#include <Eigen/Dense>

namespace input {
void gatherData(std::string,
        int &,
        std::vector<int> **,
                std::vector<double> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **,
                std::vector<std::vector<double>> **

);

void gatherData(std::string,
        int &,
        std::vector<int> **,
                std::vector<double> **,
                Eigen::MatrixXd **,
                Eigen::MatrixXd **,
                Eigen::MatrixXd **,
                Eigen::MatrixXd **,
                Eigen::MatrixXd **

);

void readVector(std::string, std::vector<std::vector<double>> **);
void readVector(std::string, std::vector<double> **, int *);
void readVector(std::string, Eigen::MatrixXd **);
void readERI(std::string, std::vector<double> **);

void readGeometry(std::string, int &, std::vector<int> **,
                  std::vector<std::vector<double>> **);
void readGeometry(std::string, int &, std::vector<int> **,
                  Eigen::MatrixXd **);

void numAtoms(std::string, int &);
void printVector(std::vector<std::vector<double>> *);
void printVector(std::vector<double> *);
void printElements(std::vector<int> *);
}

#endif
