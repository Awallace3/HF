#ifndef INPUT_HPP
#define INPUT_HPP 1

#include "stdio.h"
#include <string>
#include <vector>

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

void readVector(std::string, std::vector<std::vector<double>> **);
void readVector(std::string, std::vector<double> **);

void readGeometry(std::string, int &, std::vector<int> **,
                  std::vector<std::vector<double>> **);

void numAtoms(std::string, int &);
void printVector(std::vector<std::vector<double>> *);
void printVector(std::vector<double> *);
void printElements(std::vector<int> *);
}

#endif
