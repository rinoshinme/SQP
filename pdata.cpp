#include "pdata.h"

/*
 * A simple example
 * 2-agents model, cross road
 * agent 1: (0, 1) --> (2, 1) @ t = 1
 * agetn 2: (1, 0) --> (1, 2) @ t = 2
 * speed range: both 0-5
 *
 * The set of functions used in the calculation
 * cost function f(x)
 * equality constraint g(x)
 * separation constraint h1(x)
 * velocity constraint h2(x)
 * 
 * state variable x is of length 3 * N
 * x[ 0 -  N-1] = x1
 * x[ N - 2N-1] = x21
 * x[2N - 3N-1] = x22
 * 
 */

// # of time steps in each unified time (0-1)
static const int N = 20;
static double dt = 1.0/N;
static double sep = 0.3;
double vmin = 0.0;
double vmax = 2.0;

// x.size() == 3 * N
double fx(const Vector& x)
{
	double cost = 0;
	for (int i = 0; i < N-1; ++i)
	{
		double vel = x.GetValue(i+1) - x.GetValue(i);
		cost += vel * vel * dt;
	}
	for (int i = N; i < 3*N-1; ++i)
	{
		double vel = x.GetValue(i+1) - x.GetValue(i);
		cost += vel * vel * dt;
	}
	return cost;
} 

Vector fxPrime(const Vector& x)
{
	Vector fxp(3*N);
	for (int i = 0; i < N; ++i)
	{
		if (i == 0)
			fxp.SetValue(i, 2 * dt * (x.GetValue(i) - x.GetValue(i+1)));
		else if (i ==  N - 1)
			fxp.SetValue(i, 2 * dt * (x.GetValue(i) - x.GetValue(i-1)));
		else
			fxp.SetValue(i, 2 * dt * (2 * x.GetValue(i) - x.GetValue(i+1) - x.GetValue(i-1)));
	}

	for (int i = N; i < 3*N; ++i)
	{
		if (i == N)
			fxp.SetValue(i, 2 * dt * (x.GetValue(i) - x.GetValue(i+1)));
		else if (i == 3 * N - 1)
			fxp.SetValue(i, 2 * dt * (x.GetValue(i) - x.GetValue(i-1)));
		else
			fxp.SetValue(i, 2 * dt * (2 * x.GetValue(i) - x.GetValue(i+1) - x.GetValue(i-1)));
	}
	return fxp;
}

Vector gx(const Vector& x)
{
	// only end position equality
	double g1 = x.GetValue(N-1) - 2;
	double g2 = x.GetValue(3*N-1) - 2;
	Vector g(2);
	g.SetValue(0, g1);
	g.SetValue(1, g2);
	return g;
}

// separation constraints
Vector h1x(const Vector& x, double sep)
{
	// only N pairs 
	Vector h1x(N);
	for (int i = 0; i < N; ++i)
	{
		double x1 = x.GetValue(i);
		double x2 = x.GetValue(i+N);
		double d2 = sep * sep - (x1-1)*(x1-1) + (1-x2)*(1-x2); // <=0
		h1x.SetValue(i, d2);
	}
	return h1x;
}

// velocity constraints
Vector h2x(const Vector& x, double vmin, double vmax)
{
	Vector h2x(2*(N-1) + 2*(2*N-1));
	int idx = 0;
	for (int i = 0; i< N-1; ++i)
	{
		double h1 = x.GetValue(i) - x.GetValue(i+1) + vmin * dt;
		double h2 = x.GetValue(i+1) - x.GetValue(i) - vmax * dt;
		h2x.SetValue(idx++, h1);
		h2x.SetValue(idx++, h2);
	}

	for (int i = N; i < 3*N-1; ++i)
	{
		double h1 = x.GetValue(i) - x.GetValue(i+1) + vmin * dt;
		double h2 = x.GetValue(i+1) - x.GetValue(i) - vmax * dt;
		h2x.SetValue(idx++, h1);
		h2x.SetValue(idx++, h2);
	}
	return h2x;
}

Vector hx(const Vector& x)
{
	Vector upperhx = h1x(x, sep);
	Vector lowerhx = h2x(x, vmin, vmax);
	int n1 = upperhx.GetSize();
	int n2 = lowerhx.GetSize();
	Vector h(n1 + n2);
	for (int i = 0; i < n1; ++i)
		h.SetValue(i, upperhx.GetValue(i));
	for (int i = 0; i < n2; ++i)
		h.SetValue(i+n1, lowerhx.GetValue(i));
	return h;
}

// derivatives (Matrix)
// derivative of g
Matrix gxPrime(const Vector& x)
{
	Matrix gprime(2, 3*N);
	gprime.SetValue(0, N-1, 1);
	gprime.SetValue(1, 3*N-1, 1);
	return gprime;
}

// derivative of h1x
Matrix h1xPrime(const Vector& x)
{
	Matrix h1xprime(N, 3*N);
	for (int i = 0; i < N; ++i)
	{
		h1xprime.SetValue(i, i, -(x.GetValue(i) - 1));
		h1xprime.SetValue(i, i+N, -(x.GetValue(i+N) - 1));
	}
	return h1xprime;
}


// derivative of h2x
Matrix h2xPrime(const Vector& x)
{
	Matrix h2xprime(2*(N-1) + 2*(2*N-1), 3*N);
	int idx = 0;
	for (int i = 0; i < N-1; ++i)
	{
		h2xprime.SetValue(idx, i, 1);
		h2xprime.SetValue(idx, i+1, -1);
		idx++;
		h2xprime.SetValue(idx, i, -1);
		h2xprime.SetValue(idx, i+1, 1);
		idx++;
	}

	for (int i = N; i < 3*N-1; ++i)
	{
		h2xprime.SetValue(idx, i, 1);
		h2xprime.SetValue(idx, i+1, -1);
		idx++;
		h2xprime.SetValue(idx, i, -1);
		h2xprime.SetValue(idx, i+1, 1);
		idx++;
	}
	return h2xprime;
}

Matrix hxPrime(const Vector& x)
{
	Matrix upperhxPrime = h1xPrime(x);
	Matrix lowerhxPrime = h2xPrime(x);
	return VStack(upperhxPrime, lowerhxPrime);
	/*
	int col = upperhxPrime.GetNumCols();
	int row1 = upperhxPrime.GetNumRows();
	int row2 = lowerhxPrime.GetNumRows();

	Matrix hxp(row1+row2, col);
	for (int i = 0; i < row1; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			hxp.SetValue(i, j, upperhxPrime.GetValue(i, j));
		}
	}
	for (int i = 0; i < row2; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			hxp.SetValue(i+row1, j, lowerhxPrime.GetValue(i, j));
		}
	}
	return hxp;
	*/
}

// ========================================================================
double phix(const Vector& x, double sigma)
{
	double f = fx(x);
	Vector h = hx(x);
	Vector g = gx(x);
	Vector gn(g.GetSize());
	for (int i = 0; i < g.GetSize(); ++i)
		gn.SetValue(i, (-g.GetValue(i) > 0) ? -g.GetValue(i) : 0);
	int l0 = h.GetSize();
	int m0 = g.GetSize();
	double p;
	if (l0 == 0)
		p = f + 1.0 / sigma * gn.Norm1();
	else if (m0 == 0)
		p = f + 1.0 / sigma * h.Norm1();
	else 
		p = f + 1.0 / sigma * (h.Norm1() + gn.Norm1());
	return p;
}

double dphix(const Vector& x, double sigma, const Vector& d)
{
	Vector df = fxPrime(x);
	Vector h = hx(x);
	Vector g = gx(x);
	Vector gn(g.GetSize());
	for (int i = 0; i < g.GetSize(); ++i)
		gn.SetValue(i, (-g.GetValue(i) > 0) ? -g.GetValue(i) : 0);
	int l0 = h.GetSize();
	int m0 = g.GetSize();

	double dp;
	if (l0 == 0)
		dp = Dot(df, d) - 1.0 / sigma * gn.Norm1();
	else if (m0 == 0)
		dp = Dot(df, d) - 1.0 / sigma * h.Norm1();
	else
		dp = Dot(df, d) - 1.0 / sigma * (h.Norm1() + gn.Norm1());
	return dp;
}

Vector dlax(const Vector& x, const Vector& mu, const Vector& lam)
{
	Vector df = fxPrime(x);
	Matrix Ae = hxPrime(x);
	Matrix Ai = gxPrime(x);
	Vector dl = df - MxV(Ae, mu) - MxV(Ai, lam);
	return dl;
}

// ========================================================================
