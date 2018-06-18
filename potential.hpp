#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "molecule.hpp"

typedef struct Potential {
	double potential;
	double dmax;
	Point force;
} Potential;

Potential potential(Molecule* molecule, Point point);

#endif
