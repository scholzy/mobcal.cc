#include <cmath>

#include "constants.hpp"
#include "molecule.hpp"
#include "potential.hpp"

Potential potential(Molecule* molecule, Point point)
{
	double lj_pot = 0.0;
	Point lj_der = { 0.0, 0.0, 0.0 }, ion_dip = { 0.0, 0.0, 0.0 };
	double sum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double dmax = 0.0;

	for (auto atom: *molecule) {
		Point distance_ = { std::fabs(atom.x - point.x),
				    std::fabs(atom.y - point.y),
				    std::fabs(atom.z - point.z) };
		double distance = std::sqrt(distance_.x * distance_.x +
					    distance_.y * distance_.y +
					    distance_.z * distance_.z);
		if (distance > dmax) {
			dmax = distance;
		}

		// LJ 12-6 potential
		lj_pot += 4 * lj_well(atom) *
			((std::pow(lj_radius(atom), 12) / (std::pow(distance, 12))) -
			 (std::pow(lj_radius(atom),  6) / (std::pow(distance,  6))));

		// LJ 12-6 derivative
		double lj_der_ = 4 * lj_well(atom) *
			(( 6 * std::pow(lj_radius(atom),  6) / (std::pow(distance,  8))) -
			 (12 * std::pow(lj_radius(atom), 12) / (std::pow(distance, 14))));
		lj_der.x += lj_der_ * distance_.x;
		lj_der.y += lj_der_ * distance_.y;
		lj_der.z += lj_der_ * distance_.z;

		// Ion-dipole stuff
		double rxyz3i = atom.q / std::pow(distance, 3);
		double rxyz5i = -3.0 * atom.q / pow(distance, 5);
		ion_dip.x += rxyz3i * distance_.x;
		ion_dip.y += rxyz3i * distance_.y;
		ion_dip.z += rxyz3i * distance_.z;

		sum[0] += rxyz3i + pow(distance_.x, 2) * rxyz5i;
		sum[1] += distance_.x * distance_.y * rxyz5i;
		sum[2] += distance_.x * distance_.z * rxyz5i;
		sum[3] += rxyz3i + pow(distance_.y, 2) * rxyz5i;
		sum[4] += rxyz3i + pow(distance_.z, 2) * rxyz5i;
		sum[5] += distance_.y * distance_.z * rxyz5i;
	}

	double pot_ = lj_pot - (DIPOL * (std::pow(ion_dip.x, 2) +
					 std::pow(ion_dip.y, 2) +
					 std::pow(ion_dip.z, 2)));

	Point force = { lj_der.x - (DIPOL * ((2.0 * ion_dip.x * sum[0]) +
					     (2.0 * ion_dip.y * sum[1]) +
					     (2.0 * ion_dip.z * sum[2]))),
			lj_der.y - (DIPOL * ((2.0 * ion_dip.x * sum[1]) +
					     (2.0 * ion_dip.y * sum[3]) +
					     (2.0 * ion_dip.z * sum[4]))),
			lj_der.z - (DIPOL * ((2.0 * ion_dip.x * sum[2]) +
					     (2.0 * ion_dip.y * sum[4]) +
					     (2.0 * ion_dip.z * sum[5]))) };
	Potential pot = { pot_,
			  dmax,
			  force };
	return pot;
}
