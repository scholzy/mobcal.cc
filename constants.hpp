#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cmath>

#define HELIUM 1
#define NITROGEN 2

#define BUFFER_GAS HELIUM

static const double XE = 1.60217733e-19;
static const double XN = 6.0221367e23;
static const double XK = 1.380658e-23;
static const double XEO = 8.854187817e-12;
static const double DIPOL = 0.204956e-30 / (2.0 * 4.0 * M_PI * XEO) * XE * XE;
static const double T = 298.0;
static const double EO = 1.34e-3 * XE;
static const double RO = 3.043e-10;
static const double ROMAX = 3.043e-10;
static const double TST = XK * T / EO;

static const int IPR = 100;
/* static const int ITN = 10; */
static const int ITN = 50;

/* static const int INP = 40; */
static const int INP = 40;
/* static const int IMP = 25; */
static const int IMP = 50;
static const double CMIN = 0.0005;

static const double SW1 = 0.00005;
static const double SW2 = 0.005;
static const double DTSF1 = 0.5;
static const double DTSF2 = 0.1;

static const int INWR = 1;
static const int IFAIL = 100;
static const int IFAILC = 0;

#endif
