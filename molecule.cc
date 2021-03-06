#include "molecule.hpp"

void align_com(Molecule* molecule)
{
    Point com = { 0.0, 0.0, 0.0 };
    double mass = 0.0;

    for (auto atom : *molecule) {
        double m = 0.0;
        if (atom.m == 1.0) {
            m = 1.008;
        } else if (atom.m == 12.0) {
            m = 12.01;
        } else if (atom.m == 14.0) {
            m = 14.01;
        } else if (atom.m == 16.0) {
            m = 16.00;
        }
        com.x += atom.x * m;
        com.y += atom.y * m;
        com.z += atom.z * m;
        mass += m;
    }

    com.x /= mass;
    com.y /= mass;
    com.z /= mass;

    std::cout << "COM coords: " << com.x << "\t" << com.y << "\t" << com.z << std::endl;

    for (auto& atom : *molecule) {
        atom.x = (atom.x - com.x) * 1e-10;
        atom.y = (atom.y - com.y) * 1e-10;
        atom.z = (atom.z - com.z) * 1e-10;
        std::cout << atom.m << "\t" << atom.x << "\t" << atom.y << "\t" << atom.z << std::endl;
    }
}

Molecule rotate(Molecule molecule, double theta, double phi, double gamma)
{
    Molecule new_molecule(molecule);

    for (unsigned i = 0; i < molecule.size(); ++i) {
        Atom atom = molecule[i];
        double distance_xy = std::sqrt(atom.x * atom.x + atom.y * atom.y);
        double n_theta = theta;
        if (distance_xy != 0) {
            double o_theta = std::acos(atom.x / distance_xy);
            if (atom.y < 0.0) {
                o_theta = 2.0 * M_PI - o_theta;
            }
            n_theta += o_theta;
        }
        new_molecule[i].x = std::cos(n_theta) * distance_xy;
        new_molecule[i].y = std::sin(n_theta) * distance_xy;
    }

    for (unsigned i = 0; i < molecule.size(); ++i) {
        Atom atom = molecule[i];
        Atom new_atom = new_molecule[i];
        double distance_zy = std::sqrt(atom.z * atom.z + new_atom.y * new_atom.y);
        double n_phi = phi;
        if (distance_zy != 0) {
            double o_phi = std::acos(atom.z / distance_zy);
            if (new_atom.y < 0.0) {
                o_phi = 2.0 * M_PI - o_phi;
            }
            n_phi += o_phi;
        }
        new_molecule[i].z = std::cos(n_phi) * distance_zy;
        new_molecule[i].y = std::sin(n_phi) * distance_zy;
    }

    for (unsigned i = 0; i < molecule.size(); ++i) {
        Atom new_atom = new_molecule[i];
        double distance_xy = std::sqrt(new_atom.x * new_atom.x + new_atom.y * new_atom.y);
        double n_gamma = gamma;
        if (distance_xy != 0) {
            double o_gamma = std::acos(new_atom.x / distance_xy);
            if (new_atom.y < 0.0) {
                o_gamma = 2.0 * M_PI - o_gamma;
            }
            n_gamma += o_gamma;
        }
        new_molecule[i].x = std::cos(n_gamma) * distance_xy;
        new_molecule[i].y = std::sin(n_gamma) * distance_xy;
    }

    return new_molecule;
}

double align_x(Molecule* molecule)
{
    double distance_max = 0.0;
    int index_max = 0;

    for (unsigned i = 0; i < molecule->size(); ++i) {
        Atom atom = (*molecule)[i];
        double distance = std::sqrt(atom.x * atom.x + atom.y * atom.y + atom.z * atom.z);
        if (distance > distance_max) {
            distance_max = distance;
            index_max = i;
        }
    }

    Atom atom_max = (*molecule)[index_max];
    double distance_zy = std::sqrt(atom_max.z * atom_max.z + atom_max.y * atom_max.y);
    double phi = std::acos(atom_max.z / distance_zy) + (M_PI / 2.0);

    if (atom_max.y < 0.0) {
        phi = 2.0 * M_PI - phi;
    }

    phi = 2.0 * M_PI - phi;
    double theta = 0.0, gamma = 0.0;

    Molecule new_molecule = rotate(*molecule, theta, phi, gamma);

    Atom new_max = new_molecule[index_max];

    double distance_xy = std::sqrt(new_max.x * new_max.x + new_max.y * new_max.y);
    gamma = std::acos(new_max.x / distance_xy);
    if (new_max.y < 0.0) {
        gamma = 2.0 * M_PI - gamma;
    }
    gamma = 2.0 * M_PI - gamma;
    new_molecule = rotate(*molecule, theta, phi, gamma);

    std::cout << "Old coordinates:" << std::endl;
    for (auto atom : *molecule) {
        std::cout << atom.x << "\t"
                  << atom.y << "\t"
                  << atom.z << std::endl;
    }

    std::cout << "New coordinates:" << std::endl;
    for (auto atom : new_molecule) {
        std::cout << atom.x << "\t"
                  << atom.y << "\t"
                  << atom.z << std::endl;
    }

    for (unsigned i = 0; i < molecule->size(); ++i) {
        (*molecule)[i].x = new_molecule[i].x;
        (*molecule)[i].y = new_molecule[i].y;
        (*molecule)[i].z = new_molecule[i].z;
    }

    return distance_max;
}

double lj_well(Atom atom)
{
    double well = 0.0;

#if HELIUM
    if (atom.m == 1.0) {
        well = 0.65e-3 * XE;
    } else if (atom.m == 12.0) {
        well = 1.34e-3 * XE;
    } else if (atom.m == 14.0) {
        well = 1.34e-3 * XE;
    } else if (atom.m == 16.0) {
        well = 1.34e-3 * XE;
    }

#elif NITROGEN
    double eogas = 0.06900;
    double conve = 4.2 * 0.01036427;
    if (atom.m == 1.0) {
        well = std::sqrt(eogas * 0.0189) * conve * XE;
    } else if (atom.m == 12.0) {
        well = std::sqrt(eogas * 0.0977) * conve * XE;
    } else if (atom.m == 14.0) {
        well = std::sqrt(eogas * 0.0828) * conve * XE;
    } else if (atom.m == 16.0) {
        well = std::sqrt(eogas * 0.0558) * conve * XE;
    }

#endif

    return well;
}

double lj_radius(Atom atom)
{
    double radius = 0.0;

#if HELIUM
    if (atom.m == 1.0) {
        radius = 2.38e-10;
    } else if (atom.m == 12.0) {
        radius = 3.043e-10;
    } else if (atom.m == 14.0) {
        radius = 3.043e-10;
    } else if (atom.m == 16.0) {
        radius = 3.043e-10;
    }

#elif NITROGEN
    double rogas = 3.66;
    double convr = 0.890898718;
    if (atom.m == 1.0) {
        radius = std::sqrt(rogas * 1.2409) * convr * 1.0e-10;
    } else if (atom.m == 12.0) {
        radius = std::sqrt(rogas * 3.5814) * convr * 1.0e-10;
    } else if (atom.m == 14.0) {
        radius = std::sqrt(rogas * 4.3920) * convr * 1.0e-10;
    } else if (atom.m == 16.0) {
        radius = std::sqrt(rogas * 3.2550) * convr * 1.0e-10;
    }

#endif

    return radius;
}

double mu(Molecule* molecule)
{
    int mass = 0;
    for (auto atom : *molecule) {
        mass += atom.m;
    }
#if HELIUM
    return ((4.0026 * (double)mass) / (4.0 + (double)mass)) / (XN * 1e3);
#elif NITROGEN
    return ((28.0 * (double)mass) / (28.0 + (double)mass)) / (XN * 1e3);
#endif
}

Molecule rantate(Molecule* molecule, double* rands)
{
    double theta = rands[0] * 2.0 * M_PI;
    double phi = std::asin(rands[1] * 2.0 - 1.0) + M_PI / 2.0;
    double gamma = rands[2] * 2.0 * M_PI;

    return rotate(*molecule, theta, phi, gamma);
}

double romax(Molecule* molecule)
{
    double ro = 0.0;
    for (auto atom: *molecule) {
        double r = lj_radius(atom);
        if (r > ro) {
            ro = r;
        }
    }
    
#if NITROGEN
    ro += (1.1055e-10 / 2.0);
#endif

    return ro;
}
