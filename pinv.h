#ifndef _PINV_H
#define _PINV_H

#include "utils.h"
#include <Eigen/Dense>

Vector LeastSquare(const Matrix& matrix, const Vector& b);

#endif
