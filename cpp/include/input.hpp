#ifndef INPUT_HPP
#define INPUT_HPP 1

#include "stdio.h"
#include <string>
#include <vector>

namespace input {
void showHello();
/* void readGeometry(std::string filename, std::vector<std::string> &atom_names,
 */
/*                   std::vector<double> &x_coords, std::vector<double>
 * &y_coords, */
/*                   std::vector<double> &z_coords); */
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

void read2DVector(std::string, std::vector<std::vector<double>> **);

void readGeometry(std::string, int &, std::vector<int> **,
                  std::vector<std::vector<double>> **);

void numAtoms(std::string, int &);
void print2dVector(std::vector<std::vector<double>> *);
void printElements(std::vector<int> *);
} // namespace input

#endif
