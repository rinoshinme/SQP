#ifndef _EAMXPLE_1_H
#define _EXAMPLE_1_H

#include "utils.h"

/*
 * Problem data
 */

double fx(const Vector& x);

Vector gx(const Vector& x); // equality constraints

Vector hx(const Vector& x); // inequality constraints

Vector fxPrime(const Vector& x);

Matrix gxPrime(const Vector& x);

Matrix hxPrime(const Vector& x);


// =========================================================================
double phix(const Vector& x, double sigma);

double dphix(const Vector& x, double sigma, const Vector& d);

Vector dlax(const Vector& x, const Vector& mu, const Vector& lam);
// =========================================================================
#endif
