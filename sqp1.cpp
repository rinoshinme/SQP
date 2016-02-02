#if 1
#include <GL\glut.h>
#include <iostream>
#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "pdata.h"
#include "functions.h"

#define _SQP_DEBUG 1

#define max(a, b) (a)>(b)?(a):(b)

/* 
 * Solve the SQP problem with functions defined in example1.h
 */
void Sqpm(Vector x0, 
	      Vector mu0, 
		  Vector lam0,
		  Vector* x_result, 
		  Vector* mu_result, 
		  Vector* lam_result)
{
	int maxk = 1000;
	int n = x0.GetSize();     // 20
	int l = mu0.GetSize();    // 2
	int m = lam0.GetSize();   // 116

	double rho = 0.5;
	double eta = 0.1;
	Matrix B0(n, n);
	for (int i = 0; i < n; ++i)
		B0.SetValue(i, i, 1.0);

	Vector x(x0);
	Vector mu(mu0);
	Vector lam(lam0);
	Matrix Bk(B0);
	double sigma = 0.8;
	double epsilon1 = 1.0e-5;
	double epsilon2 = 1.0e-5;
	
	Vector hk = hx(x);
	Vector gk = gx(x);
	Vector dfk = fxPrime(x);
	Matrix Ae = hxPrime(x);
	Matrix Ai = gxPrime(x);
	Matrix Ak = VStack(Ae, Ai); // (height: l + m)
	
	int k = 0;
	Vector dk;
	while (k < maxk)
	{
#ifdef _SQP_DEBUG
		std::cout << "point 1: before subproblem starts" << std::endl;
#endif

		QPSubProblem(dfk, Bk, Ae, hk, Ai, gk, &dk, &mu, &lam);

#ifdef _SQP_DEBUG
		std::cout << "point 2: after subproblem finished" << std::endl;
#endif
		// calculate max(-gx, 0)
		Vector mgx(gk.GetSize());
		for (int i = 0; i < gk.GetSize(); ++i)
			mgx.SetValue(i, max(-gk.GetValue(i), 0));
		
		double mp1 = hk.Norm1() + mgx.Norm1();
		if (dk.Norm1() < epsilon1 && mp1 < epsilon2)
			break;

		double deta = 0.05;
		double tau = max(mu.NormInf(), lam.NormInf());
		if (sigma * (tau + deta) < 1)
			sigma = sigma;
		else
			sigma = 1.0 / (tau + 2 * deta);

		int im = 0; // Armijo search
		int mk;

#ifdef _SQP_DEBUG
		std::cout << "point 3: before Armijo search" << std::endl;
#endif
		while (im <= 20)
		{
			double rhoim = pow(rho, im);
			double p1 = phix(x + rhoim * dk, sigma);
			double p2 = phix(x, sigma);
			double p3 = eta * rhoim * dphix(x, sigma, dk);
			if (p1 - p2 < p3)
			{
				mk = im;
				break;
			}
			im++;
			if (im == 20)
				mk = 10;
		}

#ifdef _SQP_DEBUG
		std::cout << "point 4: after Armijo search" << std::endl;
#endif
		double alpha = pow(rho, mk);
		Vector x1 = x + alpha * dk;
		// update constraints
		hk = hx(x1);
		gk = gx(x1);
		dfk = fxPrime(x1);
		Ae = hxPrime(x1);
		Ai = gxPrime(x1);
		Ak = VStack(Ae, Ai);
#if 0
		Matrix pinvAk = pinv(Ak); // pinv...
		Vector lamu = MxV(pinvAk.Transpose(), dfk);
#else
		Vector lamu = LeastSquare(Ak, dfk);
#endif
		// update mu and lam
		for (int i = 0; i < l; ++i)
			mu.SetValue(i, lamu.GetValue(i));
		for (int i = 0; i < m; ++i)
			lam.SetValue(i, lamu.GetValue(i+l));
		// update Bk
		Vector sk = alpha * dk;
		Vector yk = dlax(x1, mu, lam) - dlax(x, mu, lam);

		double theta;
		if (Dot(sk, yk) > 0.2 * vBv(sk, Bk))
			theta = 1.0;
		else
			theta = 0.8 * vBv(sk, Bk) / (vBv(sk, Bk) - Dot(sk, yk));

		Vector zk = theta * yk + (1 - theta) * MxV(Bk, sk);
		Bk = Bk + Out(zk, zk) / Dot(sk, zk) - Out(MxV(Bk, sk), MxV(Bk, sk)) / vBv(sk, Bk);
		x = x1;
		k++;

#ifdef _SQP_DEBUG
		std::cout << "point 5: after update complete" << std::endl;
#endif

#if 1
		DisplayVector(x);

#endif

	}
	*x_result = x;
	*mu_result = mu;
	*lam_result = lam;
}

//=================================================================================
// main procedure

double width = 600;
double height = 600;

static Vector x0;
static Vector mu0;
static Vector lam0;
static Vector x;
static Vector mu;
static Vector lam;

#define N 20

void SolveProblem()
{
	// initialize x0, mu0 and lam0
	x0 = Vector(3 * N);
	for (int i = 0; i < 3*N; ++i)
		x0.SetValue(i, i * 0.01);
	mu0 = Vector(2);
	lam0 = Vector(6 * N - 4); // 116;
	Sqpm(x0, mu0, lam0, &x, &mu, &lam);
}

void init()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, width, 0.0, height);
	SolveProblem();
}

void display()
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5);
	glBegin(GL_POINTS);
	{
		glVertex2f((100), (100));
	}
	glEnd();
	glFlush();
}

void main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(600, 600);
	glutCreateWindow("route planning");

	init();
	glutDisplayFunc(display);
	glutMainLoop();
}
#endif
