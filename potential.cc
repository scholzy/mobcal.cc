#include <string>

#include "constants.hpp"
#include "potential.hpp"

Potential potential(Molecule* molecule, Point point)
{

#if HELIUM
    double lj_pot = 0.0;
    Point lj_der = { 0.0, 0.0, 0.0 }, ion_dip = { 0.0, 0.0, 0.0 };
    double sum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double dmax = 0.0;

    for (auto atom : *molecule) {
        Point distance_ = { std::fabs(atom.x - point.x),
            std::fabs(atom.y - point.y),
            std::fabs(atom.z - point.z) };
        double distance = std::sqrt(distance_.x * distance_.x + distance_.y * distance_.y + distance_.z * distance_.z);
        if (distance > dmax) {
            dmax = distance;
        }

        // LJ 12-6 potential
        lj_pot += 4 * lj_well(atom) * ((std::pow(lj_radius(atom), 12) / (std::pow(distance, 12))) - (std::pow(lj_radius(atom), 6) / (std::pow(distance, 6))));

        // LJ 12-6 derivative
        double lj_der_ = 4 * lj_well(atom) * ((6 * std::pow(lj_radius(atom), 6) / (std::pow(distance, 8))) - (12 * std::pow(lj_radius(atom), 12) / (std::pow(distance, 14))));
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
        sum[4] += distance_.y * distance_.z * rxyz5i;
        sum[5] += rxyz3i + pow(distance_.z, 2) * rxyz5i;
    }

    double pot_ = lj_pot - (DIPOL * (std::pow(ion_dip.x, 2) + std::pow(ion_dip.y, 2) + std::pow(ion_dip.z, 2)));
    std::cout << lj_pot / XE << std::endl;

    Point force = { lj_der.x - (DIPOL * ((2.0 * ion_dip.x * sum[0]) + (2.0 * ion_dip.y * sum[1]) + (2.0 * ion_dip.z * sum[2]))),
        lj_der.y - (DIPOL * ((2.0 * ion_dip.x * sum[1]) + (2.0 * ion_dip.y * sum[3]) + (2.0 * ion_dip.z * sum[4]))),
        lj_der.z - (DIPOL * ((2.0 * ion_dip.x * sum[2]) + (2.0 * ion_dip.y * sum[4]) + (2.0 * ion_dip.z * sum[5]))) };
    Potential pot = { pot_,
        dmax,
        force };
    return pot;

#elif NITROGEN

    double lj_pot;
    Point lj_der, ion_dip;
    double sum[6];

    double pottry[2][3], pot_mol[3];
    double dpotxtry[2][3], dpotx_mol[3], dpotytry[2][3], dpoty_mol[3], dpotztry[2][3], dpotz_mol[3];
    double qpol = 0.0;
    Point dpol;
    Point dqpol;

    double dmax = 2.0 * romax(molecule);
    double bond = 1.0976e-10;
    double Ptfn = 0.0;
    double xkT = 500.0 * XK;
    double pc = -0.4825e0;
    double pc_center = -1.0 * pc;
    double dipolxx = (1.710e-30 / (2.0 * 4.0 * M_PI * XEO)) * XE * XE;
    double dipolzz = (1.710e-30 / (2.0 * 4.0 * M_PI * XEO)) * XE * XE;
    double pot_min = 1.0e8;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            // Clear out working arrays.
            for (int k = 0; k < 6; ++k) {
                sum[k] = 0.0;
            }
            lj_pot = 0.0;
            lj_der.x = 0.0, lj_der.y = 0.0, lj_der.z = 0.0;
            ion_dip.x = 0.0, ion_dip.y = 0.0, ion_dip.z = 0.0;
            qpol = 0.0;
            dqpol.x = 0.0, dqpol.y = 0.0, dqpol.z = 0.0;

            double xc = 0.0, yc = 0.0, zc = 0.0;
            for (auto atom : *molecule) {
                if (i == 0) {
                    xc = (bond / 2.0) * (2.0 * (double)(j + 1) - 3.0);
                    dpol.x = dipolzz;
                    dpol.y = dipolxx;
                    dpol.z = dipolxx;
                } else if (i == 1) {
                    yc = (bond / 2.0) * (2.0 * (double)(j + 1) - 3.0);
                    dpol.x = dipolxx;
                    dpol.y = dipolzz;
                    dpol.z = dipolxx;
                } else if (i == 2) {
                    zc = (bond / 2.0) * (2.0 * (double)(j + 1) - 3.0);
                    dpol.x = dipolxx;
                    dpol.y = dipolxx;
                    dpol.z = dipolzz;
                }

                Point center_ = { point.x - atom.x,
                                  point.y - atom.y,
                                  point.z - atom.z };
                double center = std::sqrt(center_.x * center_.x +
                                          center_.y * center_.y +
                                          center_.z * center_.z);
                Point distance_ = { center_.x + xc,
                                    center_.y + yc,
                                    center_.z + zc };
                double distance = std::sqrt(distance_.x * distance_.x +
                                            distance_.y * distance_.y +
                                            distance_.z * distance_.z);
                if (distance < dmax) {
                    dmax = distance;
                }

                lj_pot += 4 * lj_well(atom) *
                    ((std::pow(lj_radius(atom), 12) / (std::pow(distance, 12))) -
                     (std::pow(lj_radius(atom), 6) / (std::pow(distance, 6))));

                double lj_der_ = 4 * lj_well(atom) *
                    ((6 * std::pow(lj_radius(atom), 6) / (std::pow(distance, 8))) -
                     (12 * std::pow(lj_radius(atom), 12) / (std::pow(distance, 14))));

                lj_der.x += lj_der_ * distance_.x;
                lj_der.y += lj_der_ * distance_.y;
                lj_der.z += lj_der_ * distance_.z;

                // Ion-dipole stuff
                double rxyz3i = atom.q / std::pow(center, 3);
                double rxyz5i = -3.0 * atom.q / std::pow(center, 5);
                ion_dip.x += rxyz3i * center_.x;
                ion_dip.y += rxyz3i * center_.y;
                ion_dip.z += rxyz3i * center_.z;

                // Ion-dipole derivative components
                sum[0] += rxyz3i + std::pow(center_.x, 2) * rxyz5i;
                sum[1] += center_.x * center_.y * rxyz5i;
                sum[2] += center_.x * center_.z * rxyz5i;
                sum[3] += rxyz3i + std::pow(center_.y, 2) * rxyz5i;
                sum[4] += center_.y * center_.z * rxyz5i;
                sum[5] += rxyz3i + std::pow(center_.z, 2) * rxyz5i;

                // Ion-charge ""quadrupole"" potential
                double const_k = atom.q * (XE * XE) / (4.0 * M_PI * XEO);
                qpol += pc_center * const_k / center;
                qpol += pc * const_k / distance;

                // Ion-charge ""quadrupole"" derivative
                dqpol.x -= (pc_center * const_k / std::pow(center, 3)) * center_.x;
                dqpol.y -= (pc_center * const_k / std::pow(center, 3)) * center_.y;
                dqpol.z -= (pc_center * const_k / std::pow(center, 3)) * center_.z;
                dqpol.x -= (pc * const_k / std::pow(distance, 3)) * distance_.x;
                dqpol.y -= (pc * const_k / std::pow(distance, 3)) * distance_.y;
                dqpol.z -= (pc * const_k / std::pow(distance, 3)) * distance_.z;
            }

            pottry[j][i] = lj_pot - 0.5 * ((dpol.x * ion_dip.x * ion_dip.x) + (dpol.y * ion_dip.y * ion_dip.y) + (dpol.z * ion_dip.z * ion_dip.z)) + qpol;
            dpotxtry[j][i] = lj_der.x - 0.5 * ((dpol.x * 2.0 * ion_dip.x * sum[0]) + (dpol.y * 2.0 * ion_dip.y * sum[1]) + (dpol.z * 2.0 * ion_dip.z * sum[2])) + dqpol.x;
            dpotxtry[j][i] = lj_der.y - 0.5 * ((dpol.x * 2.0 * ion_dip.x * sum[1]) + (dpol.y * 2.0 * ion_dip.y * sum[3]) + (dpol.z * 2.0 * ion_dip.z * sum[4])) + dqpol.y;
            dpotxtry[j][i] = lj_der.z - 0.5 * ((dpol.x * 2.0 * ion_dip.x * sum[2]) + (dpol.y * 2.0 * ion_dip.y * sum[4]) + (dpol.z * 2.0 * ion_dip.z * sum[5])) + dqpol.z;
        }

        pot_mol[i] = pottry[0][i] + pottry[1][i];
        if (pot_mol[i] < pot_min) {
            pot_min = pot_mol[i];
        }
        dpotx_mol[i] = dpotxtry[0][i] + dpotxtry[1][i];
        dpoty_mol[i] = dpotytry[0][i] + dpotytry[1][i];
        dpotz_mol[i] = dpotztry[0][i] + dpotztry[1][i];
    }

    for (int i = 0; i < 3; ++i) {
        double temp_pot = pot_mol[i] - pot_min;
        Ptfn += std::exp((-1.0 * temp_pot) / xkT);
    }

    double pot_ = 0.0;
    Point force = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; ++i) {
        double temp_pot = pot_mol[i] - pot_min;
        double weight = std::exp((-1.0 * temp_pot) / xkT) / Ptfn;
        pot_ += weight * pot_mol[i];
        force.x += weight * dpotx_mol[i];
        force.y += weight * dpoty_mol[i];
        force.z += weight * dpotz_mol[i];
    }

    Potential pot = { pot_, dmax, force };
    return pot;

#endif

}

