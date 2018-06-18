#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "constants.hpp"

typedef struct Point {
	double x;
	double y;
	double z;
} Point;

typedef struct Atom {
	int m;
 	double q;
	double x;
	double y;
	double z;
} Atom;

typedef std::vector<Atom> Molecule;

void align_com(Molecule* molecule);

Molecule rotate(Molecule molecule, double theta, double phi, double gamma);

double align_x(Molecule* molecule);

double lj_well(Atom atom);

double lj_radius(Atom atom);

double mu(Molecule* molecule);

Molecule rantate(Molecule* molecule, double* rands);

#endif
