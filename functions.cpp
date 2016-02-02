#include "functions.h"
#include <iostream>

#define _SUBQP_DEBUG 1

void QPSubProblem(const Vector& dfk, 
				  const Matrix& Bk,
				  const Matrix& Ae,
				  const Vector& hk,
				  const Matrix& Ai,
				  const Vector& gk,
				  Vector* d_result,
				  Vector* mu_result,
				  Vector* lam_result)
{
	int n = dfk.GetLength(); // length of x
	int l = hk.GetLength();  // # of equality constraints
	int m = gk.GetLength();  // # of inequality constraints

	double gamma = 0.05;
	double epsilon = 1.0e-4;
	double rho = 0.5;
	double sigma = 0.2;
	
	double ep0 = 0.05;
	Vector mu0(l); // automatically set all elements to 0
	Vector lam0(m);// automatically set all elements to 0
	Vector d0(n);
	for (int i = 0; i < n; ++i)
		d0.SetValue(i, 1.0);
	Vector u0(1+n+l+m);
	u0.SetValue(0, ep0);

	Vector z0(1 + n + l + m);
	for (int i = 0; i < (1 + n + l + m); ++i)
	{
		if (i == 0)
			z0.SetValue(i, ep0);
		else if (i < 1 + n)
			z0.SetValue(i, d0.GetValue(i-1));
		else if (i < 1 + n + l)
			z0.SetValue(i, mu0.GetValue(i-1-n));
		else
			z0.SetValue(i, lam0.GetValue(i-1-n-l));
	}

	int k = 0; // # of iterations
	Vector z(z0);
	double ep = ep0;
	Vector d(d0);
	Vector mu(mu0);
	Vector lam(lam0);

#ifdef _SUBQP_DEBUG
	std::cout << "subPoint 1: before loop" << std::endl;
#endif

	while (k <= 10)
	{
		Vector dh = dah(ep, d, mu, lam, dfk, Bk, Ae, hk, Ai, gk);
		if (dh.Norm() < epsilon)
			break;
		Matrix A = JacobiH(ep, d, mu, lam, dfk, Bk, Ae, hk, Ai, gk);
		Vector b = beta(ep, d, mu, lam, dfk, Bk, Ae, hk, Ai, gk, gamma) * mu0 - dh;
		// inverse of A
#if 1
		Matrix invA = A.Inverse();
		Vector dz = MxV(invA, b);
#else
		Vector dz = LeastSquare(A, b); // (1 + n + l + m)
#endif
		double de;
		Vector dd(n);
		Vector du(l);
		Vector dl(m);
			
		de = dz.GetValue(0);
		for (int i = 0; i < n; ++i)
			dd.SetValue(i, dz.GetValue(i+1));
		for (int i = 0; i < l; ++i)
			du.SetValue(i, dz.GetValue(i+1+n));
		for (int i = 0; i < m; ++i)
			dl.SetValue(i, dz.GetValue(i+1+n+l));

		int i = 0;
		int mk = 0;
		int mm = 0;
		Vector dh1;

#ifdef _SUBQP_DEBUG
		std::cout << "subPoint 2: iteration begin" << std::endl;
#endif

#if 1
		while (mm <= 20)
		{
			double rhoi = pow(rho, i);
			dh1 = dah(ep+rhoi*de, d+rhoi*dd, mu+rhoi*du, lam+rhoi*dl, dfk, Bk, Ae, hk, Ai, gk);

			if (dh1.Norm() <= (1 - sigma * (1 - gamma * ep0) * rhoi) * dh.Norm())
			//if (dh1.Norm() < dh.Norm())
			{
				mk = i;
				break;
			}
			i = i + 1;
			if (i == 20)
				mk = 10;
		}
#endif

#ifdef _SUBQP_DEBUG
		std::cout << "subPoint 2: iteration end" << std::endl;
#endif
		double alpha = pow(rho, mk);

		ep = ep + alpha * de;
		d = d + alpha * dd;
		mu = mu + alpha * du;
		lam = lam + alpha * dl;
		k++;
		// std::cout << "k = " << k << std::endl;
	}
#ifdef _SUBQP_DEBUG
	std::cout << "subPoint 3: end subproblem" << std::endl;
#endif

	*d_result = d;
	*mu_result = mu;
	*lam_result = lam;
}

double phi(double ep, double a, double b)
{
	return a + b - sqrt(a * a + b * b + 2 * ep * ep);
}

Vector dah(double ep,
		   const Vector& d,
		   const Vector& mu,
		   const Vector& lam,
		   const Vector& dfk,
		   const Matrix& Bk,
		   const Matrix& Ae,
		   const Vector& hk,
		   const Matrix& Ai,
		   const Vector& gk)
{
	int n = dfk.GetSize();
	int l = hk.GetSize();
	int m = gk.GetSize();

	Vector dh(n+l+m+1); // automatically set to 0
	
	// length n
	Vector r1 = MxV(Bk, d) - MxV(Ae.Transpose(), mu) - MxV(Ai.Transpose(), lam) + dfk;
	// length l (# of equality constraints)
	Vector r2 = hk + MxV(Ae, d);

	dh.SetValue(0, ep);
	for (int i = 0; i < n; ++i)
		dh.SetValue(i + 1, r1.GetValue(i));
	for (int i = 0; i < l; ++i)
		dh.SetValue(i + 1 + n, r2.GetValue(i));
	for (int i = 0; i < m; ++i)
	{
		double sum = 0.0; // Ai(i, :) * d
		for (int k = 0; k < d.GetSize(); ++k)
			sum += Ai.GetValue(i, k) * d.GetValue(k);
		double p = phi(ep, lam.GetValue(i), gk.GetValue(i) + sum);
		dh.SetValue(i + 1 + n + l, p);
	}
	return dh;
}

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
			double gamma)
{
	Vector dh = dah(ep, d, mu, lam, dfk, Bk, Ae, hk, Ai, gk);
	double n = dh.Norm();
	double m = (1 < n) ? 1 : n;
	return gamma * n * m;
}

Vector ddv(double ep,
		   const Vector& d,
		   const Vector& lam,
		   const Matrix& Ai,
		   const Vector& gk,
		   Matrix* dd1,
		   Matrix* dd2)
{
	int m = gk.GetSize();
	dd1->Resize(m, m);
	dd2->Resize(m, m);
	Vector v1(m);

	for (int i = 0; i < m; ++i)
	{
		// calculate Ai(i, :) * d
		double sum = 0.0;
		for (int k = 0; k < d.GetSize(); ++k)
			sum += Ai.GetValue(i, k) * d.GetValue(k);

		double fm2 = lam.GetValue(i) * lam.GetValue(i) + 
			(gk.GetValue(i) + sum) * (gk.GetValue(i) + sum) + 
			2 * ep * ep;
		double fm = sqrt(fm2);
		dd1->SetValue(i, i, 1 - lam.GetValue(i) / fm);
		dd2->SetValue(i, i, 1 - (gk.GetValue(i) + sum) / fm);
		v1.SetValue(i, -2 * ep / fm);
	}
	return v1;
}

Matrix JacobiH(double ep,
			   const Vector& d,
			   const Vector& mu,
			   const Vector& lam,
			   const Vector& dfk,
			   const Matrix& Bk,
			   const Matrix& Ae,
			   const Vector& hk,
			   const Matrix& Ai,
			   const Vector& gk)
{
	int n = dfk.GetSize();
	int l = hk.GetSize();
	int m = gk.GetSize();

	Matrix A(1+n+l+m, 1+n+l+m);
	Matrix dd1;
	Matrix dd2;
	Vector v1 = ddv(ep, d, lam, Ai, gk, &dd1, &dd2);

	Matrix dd2Ai;
	MxM(dd2, Ai, &dd2Ai);
	A.SetValue(0, 0, 1.0);

	// Bk
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			A.SetValue(i + 1, j + 1, Bk.GetValue(i, j));

	// -Ae'
	Matrix Aet = Ae.Transpose();
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < l; ++j)
			A.SetValue(i + 1, j + 1 + n, -Aet.GetValue(i, j));

	// -Ai'
	Matrix Ait = Ai.Transpose();
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			A.SetValue(i + 1, j + 1 + n + l, -Ait.GetValue(i, j));

	// Ae
	for (int i = 0; i < l; ++i)
		for (int j = 0; j < n; ++j)
			A.SetValue(i + 1 + n, j + 1, Ae.GetValue(i, j));

	// v1
	for (int i = 0; i < m; ++i)
		A.SetValue(i + 1 + n + l, 0, v1.GetValue(i));

	// dd2Ai
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			A.SetValue(i + 1 + n + l, j + 1, dd2Ai.GetValue(i, j));

	// dd1
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < m; ++j)
			A.SetValue(i + 1 + n + l, i + 1 + n + l, dd1.GetValue(i, j));
	
	return A;
}
