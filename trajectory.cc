#include <string>

#include "trajectory.hpp"

void derivative(Molecule* molecule, double* w, double* dw)
{
    dw[0] = w[1] / mu(molecule);
    dw[2] = w[3] / mu(molecule);
    dw[4] = w[5] / mu(molecule);

    Point point = { w[0], w[2], w[4] };
    Potential p = potential(molecule, point);
    dw[1] = -1.0 * p.force.x;
    dw[3] = -1.0 * p.force.y;
    dw[5] = -1.0 * p.force.z;
}

void diffeq(Molecule* molecule, double* w, double* dw, double* dt, int* l, double* q, double* hvar, double* hcvar, double* time, double array[][40])
{
    double a[4] = { 0.50, 0.292893218814, 1.70710678118, 0.1666666666667 };
    double b[4] = { 2.0, 1.0, 1.0, 2.0 };
    double c[4] = { -0.5, -0.292893218814, -1.70710678118, -0.5 };

    double ampc[5] = { -0.111059153612,
        0.672667757774,
        -1.70633621697,
        2.33387888707,
        -1.8524668225 };
    double amcc[4] = { 0.0189208128941,
        -0.121233356692,
        0.337771548703,
        -0.55921513665 };
    double var = 2.97013888888;
    double cvar = 0.990972222222;
    double acst = 0.332866152768;

    if (*l == 0) {
        for (int j = 0; j < 6; ++j) {
            q[j] = 0.0;
        }
        *hvar = (*dt) * var;
        *hcvar = (*dt) * cvar;
    }

    if (*l >= 0) {
        *l += 1;
        int k = 0;
        do {
            for (int j = 0; j < 4; j++) {
                if (std::pow(-1, j + 1) > 0) {
                    *time += 0.5 * (*dt);
                }
                for (int i = 0; i < 6; i++) {
                    dw[i] *= (*dt);
                    double r = a[j] * (dw[i] - b[j] * q[i]);
                    // std::cout << r << std::endl;
                    w[i] += r;
                    q[i] += 3.0 * r + c[j] * dw[i];
                }
            }

            derivative(molecule, w, dw);
            k += 1;
        } while (k != 1);

        if ((*l) - 6 >= 0) {
            *l = -1;
            *dt *= 2.0;
            return;
        } else {
            for (int j = 0; j < 6; j++) {
                array[(*l) - 1][j] = dw[j];
            }
            return;
        }
    } else {
        double savw[6] = { w[0], w[1], w[2], w[3], w[4], w[5] };
        double savdw[6] = { dw[0], dw[1], dw[2], dw[3], dw[4], dw[5] };
        for (int j = 0; j < 6; j++) {
            array[5][j] = savdw[j];
            for (int i = 0; i < 5; i++) {
                array[5][j] += ampc[i] * array[i][j];
            }
            w[j] += array[5][j] * (*hvar);
        }
        *time += *dt;
        derivative(molecule, w, dw);
        for (int j = 0; j < 6; j++) {
            array[5][j] = acst * dw[j];
            for (int i = 0; i < 4; i++) {
                array[i][j] = array[i + 1][j];
                array[5][j] += array[i][j] * amcc[i];
            }
            array[4][j] = savdw[j];
            w[j] = savw[j] + (*hcvar) * (array[4][j] + array[5][j]);
        }
        derivative(molecule, w, dw);
    }
}

Trajout trajectory(Molecule* molecule, double v, double b)
{
    Point vel = { 0.0, -1.0 * v, 0.0 };

    /* Determine time step for integration */
    double top = (v / 95.2381) - 0.5;
    if (v >= 1000.0) {
        top = 10.0;
    } else if (v >= 2000) {
        top = 10.0 - ((v - 2000.0) * 7.5e-3);
    } else if (v >= 3000.0) {
        top = 2.5;
    }

    double dt1 = top * DTSF1 * 1.0e-11 / v;
    double dt2 = dt1 * DTSF2;

    double e0 = 0.5 * mu(molecule) * std::pow(v, 2);
    Point coord = { b, 0.0, 0.0 };

    double y_max = 0.0;
    double y_min = 0.0;
    for (auto atom : *molecule) {
        if (atom.y > y_max) {
            y_max = atom.y;
        }
        if (atom.y < y_min) {
            y_min = atom.y;
        }
    }
    y_max += 1e-10;
    y_min -= 1e-10;

    coord.y = y_max;
    Potential p = potential(molecule, coord);

    if (std::fabs(p.potential / e0) <= SW1) {
        do {
            coord.y -= 1e-10;
            p = potential(molecule, coord);
        } while (std::fabs(p.potential) / e0 < SW1);
    } else {
        do {
            coord.y += 10e-10;
            p = potential(molecule, coord);
        } while (std::fabs(p.potential) / e0 > SW1);
        do {
            coord.y -= 1e-10;
            p = potential(molecule, coord);
        } while (std::fabs(p.potential) / e0 > SW1);
    }

    double e_total = e0 + p.potential;

    /* Trajectory is ready to go! */

    double w[6] = { coord.x, vel.x * mu(molecule),
        coord.y, vel.y * mu(molecule),
        coord.z, vel.z * mu(molecule) };
    double dw[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    derivative(molecule, w, dw);

    int l = 0;
    int ns = 0;
    int nw = 0;
    double dt = dt1;
    double hvar = 0.0;
    double hcvar = 0.0;
    double q[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double time = 0.0;
    double array[6][40];

    do {
        do {
            diffeq(molecule, w, dw, &dt, &l, q, &hvar, &hcvar, &time, array);
            // std::cout << w[0] << " " << w[2] << " " << w[4] << std::endl;
            nw += 1;
        } while (nw != INWR);

        ns += nw;
        nw = 0;

        if (ns > 30000) {
            break;
        }

        Point point = { w[0],
            w[2],
            w[4] };
        p = potential(molecule, point);

        if (p.dmax < ROMAX) {
            continue;
        }

        if ((std::fabs(p.potential) / e0 > SW2) && (dt == dt1)) {
            dt = dt2;
            l = 0;
        }

        if ((std::fabs(p.potential) / e0 < SW2) && (dt == dt2)) {
            dt = dt1;
            l = 0;
        }

        if (std::fabs(p.potential) / e0 > SW1) {
            continue;
        }

        if (ns < 50) {
            continue;
        }

        break;
    } while (1);

    double num = 0.0;
    double den = 0.0;
    double ang = 0.0;

    if (dw[0] >= 0) {
        num = dw[2] * -1.0 * v;
        den = v * std::sqrt(std::pow(dw[0], 2) + std::pow(dw[2], 2) + std::pow(dw[4], 2));
        ang = std::acos(num / den);
    } else {
        num = dw[2] * -1.0 * v;
        den = v * std::sqrt(std::pow(dw[0], 2) + std::pow(dw[2], 2) + std::pow(dw[4], 2));
        ang = -1.0 * std::acos(num / den);
    }

    double e = 0.5 * mu(molecule) * (std::pow(dw[0], 2) + std::pow(dw[2], 2) + std::pow(dw[4], 2));
    Trajout output = { ang,
        (e + p.potential) / e_total };
    return output;
}
