#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "utils.h"

/*
 * find the solution of the qp subproblem 
 *     f = 0.5 * d * Bk * d + dfk * d;
 * subject to 
 *     hk + Ae * d = 0
 *     gk + Ai * d <= 0
 * 
 */
void QPSubProblem(const Vector& dfk, 
				  const Matrix& Bk,
				  const Matrix& Ae,
				  const Vector& hk,
				  const Matrix& Ai,
				  const Vector& gk,
				  Vector* d_result,
				  Vector* mu_result,
				  Vector* lam_result);

double phi(double ep, double a, double b);

Vector dah(double ep,
		 const Vector& d,
		 const Vector& mu,
		 const Vector& lam,
		 const Vector& dfk,
		 const Matrix& Bk,
		 const Matrix& Ae,
		 const Vector& hk,
		 const Matrix& Ai,
		 const Vector& gk);

double beta(double ep,
			const Vector& d,
			const Vector& mu,
			const Vector& lam,
			const Vector& dfk,
			const Matrix& Bk,
			const Matrix& Ae,
			const Vector& hk,
			const Matrix& Ai,
			const Vector& gk,
			double gamma);

Vector ddv(double ep,
		 const Vector& d,
		 const Vector& lam,
		 const Matrix& Ai,
		 const Vector& gk,
		 Matrix* dd1,
		 Matrix* dd2);


Matrix JacobiH(double ep,
			   const Vector& d,
			   const Vector& mu,
			   const Vector& lam,
			   const Vector& dfk,
			   const Matrix& Bk,
			   const Matrix& Ae,
			   const Vector& hk,
			   const Matrix& Ai,
			   const Vector& gk);

#endif
