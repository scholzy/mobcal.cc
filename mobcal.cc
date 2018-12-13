#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "constants.hpp"
#include "molecule.hpp"
#include "potential.hpp"
#include "trajectory.hpp"

Molecule get_input(std::string file)
{
    std::fstream input_file(file.c_str());
    double elem, x, y, z, charge;
    Molecule molecule;

    while (input_file >> elem >> x >> y >> z >> charge) {
        Atom a = { elem,
            charge,
            x,
            y,
            z };
        molecule.push_back(a);
    }

    return molecule;
}

Point get_extents(Molecule* molecule, double dmax)
{
    int irn = 1000;
    double ddd = (romax(molecule) + dmax) / (double)irn;
    std::cout << ddd << std::endl;
    Point e_max = { 0.0, 0.0, 0.0 };
    Point r_max = { 0.0, 0.0, 0.0 };
    Point r00_max = { 0.0, 0.0, 0.0 };
    Point coord = { 0.0, 0.0, 0.0 };

    for (int i = 1; i <= irn; ++i) {
        coord.x = romax(molecule) + dmax - ((double)i * ddd);
        coord.y = 0.0, coord.z = 0.0;
        Potential p = potential(molecule, coord);
        if (p.potential > 0.0) {
            break;
        }
        r00_max.x = coord.x;
        if (p.potential < e_max.x) {
            e_max.x = p.potential;
            r_max.x = coord.x;
        }
    }

    std::cout << "Y STARTING" << std::endl;

    for (int i = 1; i <= irn; ++i) {
        coord.y = romax(molecule) + dmax - ((double)i * ddd);
        coord.x = 0.0, coord.z = 0.0;
        Potential p = potential(molecule, coord);
        /* std::cout << coord.y/1e-10 << "\t" << p.potential/XE << std::endl; */
        if (p.potential > 0.0) {
            break;
        }
        r00_max.y = coord.y;
        if (p.potential < e_max.y) {
            e_max.y = p.potential;
            r_max.y = coord.y;
        }
    }

    for (int i = 1; i <= irn; ++i) {
        coord.z = romax(molecule) + dmax - ((double)i * ddd);
        coord.x = 0.0, coord.y = 0.0;
        Potential p = potential(molecule, coord);
        if (p.potential > 0.0) {
            break;
        }
        r00_max.z = coord.z;
        if (p.potential < e_max.z) {
            e_max.z = p.potential;
            r_max.z = coord.z;
        }
    }

    std::cout << "Well depths along the x-, y-, and z-axes:" << std::endl;
    std::cout << e_max.x / XE << "\t"
              << e_max.y / XE << "\t"
              << e_max.z / XE << "\n"
              << std::endl;

    std::cout << "Well coordinates along the x-, y-, and z-axes:" << std::endl;
    std::cout << r_max.x / 1e-10 << "\t"
              << r_max.y / 1e-10 << "\t"
              << r_max.z / 1e-10 << "\n"
              << std::endl;

    std::cout << "Distance where e = 0 along the x-, y-, and z-axes:" << std::endl;
    std::cout << r00_max.x / 1e-10 << "\t"
              << r00_max.y / 1e-10 << "\t"
              << r00_max.z / 1e-10 << "\n"
              << std::endl;

    return r_max;
}

void setup_gst(Molecule* molecule, double* wgst, double* pgst)
{
    double dgst = 5.0e-7 * 6.0 * sqrt(TST);
    double gst = dgst;
    double sum = 0.0, sum1 = 0.0, sum2 = 0.0;
    double hold1 = 0.0, hold2 = 0.0, hold3 = 0.0;

    for (int i = 0; i < INP; ++i) {
        sum1 += sqrt((double)i + 1.0);
    }

    for (int i = 0; i < INP; ++i) {
        hold1 = sqrt((double)i + 1.0);
        hold2 = sqrt((double)i);
        sum2 += hold2;
        wgst[i] = hold1 / sum1;
        double gstt = std::pow(TST, 3) * (sum2 + (hold1 / 2.0)) / sum1;

        do {
            sum += std::exp(-1.0 * std::pow(gst, 2) / TST) * std::pow(gst, 5) * dgst;
            gst += dgst;
            if (sum > gstt) {
                pgst[i] = gst - (dgst / 2.0);
            }
        } while (sum < gstt);

        hold1 = std::sqrt((pgst[i] * pgst[i] * EO) / (0.5 * mu(molecule)));
        hold2 = 0.5 * mu(molecule) * std::pow(hold1, 2) / (XK * T);
        hold3 = std::exp(-1.0 * std::pow(pgst[i], 2) / TST) * std::pow(pgst[i], 5);

        std::cout << pgst[i] << "\t" << wgst[i] << "\t" << hold1 << "\t" << hold2 << "\t" << hold3 << "\t" << sum / std::pow(TST, 3) << std::endl;
    }
}

void setup_bst(Molecule* molecule, double* wgst, double* pgst, double* b2max, double* cosx, Point extents)
{
    double dbst2 = 1.0;
    double dbst22 = dbst2 / 10.0;

    for (int i = 0; i < INP; ++i) {
        b2max[i] = 0.0;
    }

    for (int i = INP - 1; i >= 0; --i) {
        double gst2 = pgst[i] * pgst[i];
        double v = std::sqrt((gst2 * EO) / (0.5 * mu(molecule)));
        int ibst = std::trunc(extents.x / RO) - 6;

        if (i < INP - 1) {
            ibst = std::trunc(b2max[i+1] / dbst2) - 6;
        }

        if (ibst < 0) {
            ibst = 0;
        }

        double bst2, b;
        Trajout trj;
        do {
            bst2 = dbst2 * ((double)ibst);
            b = RO * std::sqrt(bst2);
            trj = trajectory(molecule, v, b);
            cosx[ibst] = 1.0 - std::cos(trj.ang);
            /* std::cout << b << "\t" << bst2 << "\t" << trj.ang << "\t" << cosx[ibst] << "\t" << trj.erat << std::endl; */

            if (ibst >= 5) {
                if ((cosx[ibst] < CMIN) && (cosx[ibst - 1] < CMIN) && (cosx[ibst - 2] < CMIN) && (cosx[ibst - 3] < CMIN) && (cosx[ibst - 4] < CMIN)) {
                    break;
                }
            }
            ibst += 1;
        } while (1);

        b2max[i] = (double)(ibst - 5) * dbst2;
        do {
            b2max[i] += dbst22;
            b = RO * std::sqrt(b2max[i]);
            trj = trajectory(molecule, v, b);
            if (1.0 - std::cos(trj.ang) < CMIN) {
                break;
            }
        } while (1);
    }

    for (int i = 0; i < INP; ++i) {
        std::cout << pgst[i] << "\t" << b2max[i] << "\t" << RO * std::sqrt(b2max[i]) * 1.0e10 << std::endl;
    }
}

void calculate(Molecule* molecule, double* wgst, double* pgst, double* b2max, double* cosx, std::mt19937_64 rng)
{
    double q1st[INP] = { 0.0 };
    double q2st[INP] = { 0.0 };

    double om11st[ITN] = { 0.0 };
    double om12st[ITN] = { 0.0 };
    double om13st[ITN] = { 0.0 };
    double om22st[ITN] = { 0.0 };

    std::uniform_real_distribution<double> unif;
    for (int i = 0; i < ITN; ++i) {
        om11st[i] = 0.0;
        om12st[i] = 0.0;
        om13st[i] = 0.0;
        om22st[i] = 0.0;

        for (int j = 0; j < INP; ++j) {
            double gst2 = std::pow(pgst[j], 2);
            double v = std::sqrt((gst2 * EO) / (0.5 * mu(molecule)));
            double temp1 = 0.0, temp2 = 0.0;

#pragma omp parallel for
            for (int k = 0; k < IMP; ++k) {
                double rnb = unif(rng);
                double rands[3] = { unif(rng), unif(rng), unif(rng) };
                Molecule mol = rantate(molecule, rands);
                double bst2 = rnb * b2max[j];
                double b = RO * std::sqrt(bst2);
                Trajout trj = trajectory(&mol, v, b);
                double hold1 = 1.0 - std::cos(trj.ang);
                double hold2 = std::pow(std::sin(trj.ang), 2);
                temp1 += hold1 * b2max[j] / (double)IMP;
                temp2 += 1.5 * hold2 * b2max[j] / (double)IMP;
            }

            om11st[i] += temp1 * wgst[j];
            om12st[i] += temp1 * pgst[j] * pgst[j] * wgst[j] * (1.0 / (3.0 * TST));
            om13st[i] += temp1 * pow(pgst[j], 4) * wgst[j] * (1.0 / (12.0 * TST * TST));
            om22st[i] += temp2 * pgst[j] * pgst[j] * wgst[j] / (1.0 / (3.0 * TST));
            q1st[j] += temp1;
            q2st[j] += temp2;
        }

        printf("Cycle %2d: %f\n", i, om11st[i] * M_PI * RO * RO * 1e20);
    }

    double sum = 0.0;
    for (int i = 0; i < ITN; i++) {
        sum += om11st[i] * M_PI * RO * RO * 1e20;
    }
    sum /= (double)ITN;
    printf("Average CCS: %f\n", sum);
}

void mobil2(Molecule molecule)
{
    const unsigned int seed = time(0);
    std::mt19937_64 rng(seed);

    align_com(&molecule);

    double distance_max = align_x(&molecule);

    Point extents = get_extents(&molecule, distance_max);

    double wgst[INP];
    double pgst[INP];
    setup_gst(&molecule, wgst, pgst);

    std::cout << "Done with gst!" << std::endl;

    double b2max[INP];
    double cosx[2000];
    setup_bst(&molecule, wgst, pgst, b2max, cosx, extents);

    std::cout << "Done with bst!" << std::endl;

    calculate(&molecule, wgst, pgst, b2max, cosx, rng);
}

int main(int argc, char** argv)
{
    auto input = (std::string)argv[1];
    Molecule molecule = get_input(input);

    mobil2(molecule);

    return 0;
}
