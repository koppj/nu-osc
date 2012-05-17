// Definitions of physical constants
#ifndef __CONST_H
#define __CONST_H

// Mathematics
#define M_SQRT_PI     1.772453850905516027298 // sqrt(pi)

// Particle physics
#define M_EL          0.51099892e6// [eV] Electron mass
#define M_MU        105.658369e6  // [eV] Muon mass
#define M_TAU      1776.99e6      // [eV] Tau mass
#define M_P         938.3e6       // [eV] Proton mass
#define M_PI_PLUS   139.57018e6   // [eV] Charged pion mass
#define M_PI_MINUS  M_PI_PLUS
#define M_PI_0      134.9766e6    // [eV] Neutral pion mass
#define M_U         931.494e6     // [eV] Atomic mass unit
#define ALPHA_0     (1/137.036)   // Fine structure constant at zero energy
#define ALPHA_MZ    (1/128.0)     // Fine structure constant at M_Z
#define ALPHA_10GEV ALPHA_MZ      // Fine structure constant at 10 GeV (approximately)
#define M_Z          91.1876e9    // [eV] mass of Z boson

// Astrophysics
#define C           2.9979e8      // [m/s] speed of light
#define V_SUN    (220.0e3/C)      // Velocity of Sun relative to GC (in units of c)
#define V_EARTH   (29.8e3/C)      // Velocity of Earth relative to Sun (in units of c)
#define V_ESC    (550.0e3/C)      // Galactic escape velocity
#define XI          1.0           // 100% of halo composed of chi0
#define RHO_0       0.3 * GEV / (CM*CM*CM) // WIMP density

// Unit conversion
#define PICOBARN  (1.0e-36 * SQR(CM))// [eV^-2 / pb]
#define FEMTOBARN (1.0e-39 * SQR(CM))// [eV^-2 / fb]

#define TEV        1.0e12            // [eV/TeV]
#define GEV        1.0e9             // [eV/GeV]
#define MEV        1.0e6             // [eV/MeV]
#define KEV        1.0e3             // [eV/keV]
#define JOULE      (1/1.60225e-19)   // [eV/J]

#define KG         5.62e35           // [eV/kg]
#define GRAMS      (1e-3*KG)         // [eV/g]

#define SEC        1.523e15          // [eV^-1 / s]
#define HOURS      (3600.*SEC)       // [eV^-1 / h]
#define DAYS       (24.*HOURS)       // [eV^-1 / d]
#define YRS        (365.*DAYS)       // [eV^-1 / yr]

#define T_JUN_2ND  (152.5/365.*YRS)  // ~ 2nd June (V_EARTH aligned with V_SUN)
#define T_DEC_2ND  (335./365.*YRS)   // ~ 2nd December (V_EARTH anti-aligned with V_SUN)
#define T_MEAN     (0.668*YRS)       // Mean between June 2nd and December 2nd
                                     // --> expected zeros of DAMA oscillations

#define KM         (1.e3 * METER)
#define METER      5.076e6           // [eV^-1 / m]
#define CM         (1.e-2 * METER)
#define FERMI      (1.0e-15*METER)   // [eV^-1 / fm]
#define ANGSTROM   (1.e-10 * METER)

#define PA         (JOULE/(METER*METER*METER)) // [eV^4/Pa]
#define HPA        (100.0*PA)        // [eV^4/hPa]
#define ATM        (101325*PA)       // [eV^4/atm]
#define PSI        (6895*PA)         // [eV^4/psi (pound per square inch)]

#define KELVIN     (1/1.1604505e4)   // [eV / K]

#define DEGREES    (M_PI/180.0)      // [rad/degree]
#define DEG        (DEGREES)         // [rad/degree]

/* Effective electron numbers for calculation of matter profile */
#define Ne_MANTLE       0.497
#define Ne_CORE         0.468

#define REARTH           6371  /* km */
#define RCORE            3480  /* km */

#endif

