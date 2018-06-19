#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include <string>

#include "constants.hpp"
#include "molecule.hpp"
#include "potential.hpp"

typedef struct Trajout {
    double ang;
    double erat;
} Trajout;

void derivative(Molecule* molecule, std::string gas, double* w, double* dw);

void diffeq(Molecule* molecule, std::string gas, double* w, double* dw, double* dt, int* l, double* q, double* hvar, double* hcvar, double* time, double array[][40]);

Trajout trajectory(Molecule* molecule, std::string gas, double v, double b);

#endif
