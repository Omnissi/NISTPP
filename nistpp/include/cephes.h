#pragma once

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                   C E P H E S   M A T H   F U N C T I O N S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

double cephes_igamc(double a, double x);
double cephes_igam(double a, double x);
double cephes_lgam(double x);
double cephes_polevl(double x, double *coef, int N);
double cephes_p1evl(double x, double *coef, int N);
double cephes_erf(double x);
double cephes_erfc(double x);
double cephes_normal(double x);
