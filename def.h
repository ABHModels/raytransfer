#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

#define imax 100 //Number of radial values

/* RK45 constants */
#define a1 1.0/4.0
#define b1 3.0/32.0
#define b2 9.0/32.0
#define c1 1932.0/2197.0
#define c2 -7200.0/2197.0
#define c3 7296.0/2197.0
#define d1 439.0/216.0
#define d2 -8.0
#define d3 3680.0/513.0
#define d4 -845.0/4104.0
#define e1 -8.0/27.0
#define e2 2.0
#define e3 -3544.0/2565.0
#define e4 1859.0/4104.0
#define e5 -11.0/40.0
#define f1 25.0/216.0
#define f2 0.0
#define f3 1408.0/2565.0
#define f4 2197.0/4104.0
#define f5 -1.0/5.0
#define g1 16.0/135.0
#define g2 0.0
#define g3 6656.0/12825.0
#define g4 28561.0/56430.0
#define g5 -9.0/50.0
#define g6 2.0/55.0

#define Pi acos(-1.0)

/* global variables to avoid passing to functions */
long double xscr, yscr;
long double defpar, epsi3, a13, a22, a52;
long double spin, inc, isco;
long double spin2 = spin*spin;
long double Mdl, eta;

void xyfromrphi(long double rscr, long double pscr, long double rdisk);
void raytrace(long double xscr, long double yscr, long double traced[], long double rdisk);
void rayprecise(long double rdisk, long double germtol, long double pscr, long double traced[]);
void diffeqs(long double b, long double vars[], long double diffs[]);
void redshift(long double r, long double th, long double ktkp, long double& gg);
long double specific_energy(long double r);
long double emis_angle(long double r, long double th, long double kr, long double kth);
void metric(long double r, long double th, long double g[][4]);
void metric_rderivatives(long double r, long double th, long double dg[][4]);
void metric_r2derivatives(long double r, long double th, long double dg2[][4]);
void uppermetric(long double r, long double th, long double rth[]);
void find_isco();
void gauleg(long double rdisk_i, long double rdisk_f, long double rdisk[]);

#include "diffeqs.cpp"
#include "rayprecise.cpp"
#include "metric.cpp"
#include "raytracing.cpp"
#include "redshift.cpp"
#include "findisco.cpp"
#include "gauleg.cpp"

#endif
