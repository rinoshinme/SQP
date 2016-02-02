#if 1
#include <GL\glut.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "pdata.h"
#include "functions.h"
#include "pinv.h"

/* 
 * Solve the SQP problem with functions defined in example1.h
 */
void Sqpm(Vector x0, Vector mu0, Vector lam0,
		  Vector* x_result, Vector* mu_result, Vector* lam_result)
{
	int maxk = 1000;
	int n = x0.GetSize();
	int l = mu0.GetSize();
	int m = lam0.GetSize();

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
	double epsilon1 = 1.0e-6;
	double epsilon2 = 1.0e-5;

	Vector hk = hx(x);
	Vector gk = gx(x);
	Vector dfk = fxPrime(x);
	Matrix Ae = hxPrime(x);
	Matrix Ai = gxPrime(x);
	Matrix Ak = VStack(Ae, Ai);

	int k = 0;
	Vector dk;
	while (k < maxk)
	{
		QPSubProblem(dfk, Bk, Ae, hk, Ai, gk, &dk, &mu, &lam);
		// calculate max(-gx, 0)
		Vector mgx(gk.GetSize());
		for (int i = 0; i < gk.GetSize(); ++i)
			mgx.SetValue(i, (-gk.GetValue(i) > 0) ? -gk.GetValue(i) : 0);
		
		double mp1 = hk.Norm1() + mgx.Norm1();
		if (dk.Norm1() < epsilon1 && mp1 < epsilon2)
			break;

		double deta = 0.05;
		double tau = (mu.NormInf() > lam.NormInf()) ? mu.NormInf() : lam.NormInf();
		if (sigma * (tau + deta) < 1)
			sigma = sigma;
		else
			sigma = 1.0 / (tau + 2 * deta);

		int im = 0; // Armijo search
		int mk;
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
		double alpha = pow(rho, mk);
		Vector x1 = x + alpha * dk;
		// update
		hk = hx(x1);
		gk = gx(x1);
		dfk = fxPrime(x1);
		Ae = hxPrime(x);
		Ai = gxPrime(x);
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
		Vector yk = dlax(x1, mu, lam) -dlax(x, mu, lam);

		double theta;
		if (Dot(sk, yk) > 0.2 * vBv(sk, Bk))
			theta = 1.0;
		else
			theta = 0.8 * vBv(sk, Bk) / (vBv(sk, Bk) - Dot(sk, yk));

		Vector zk = theta * yk + (1-theta) * MxV(Bk, sk);
		Bk = Bk + Out(zk, zk) / Dot(sk, zk) - Out(MxV(Bk, sk), MxV(Bk, sk)) / vBv(sk, Bk);
		x = x1;
		k++;
	}
	*x_result = x;
	*mu_result = mu;
	*lam_result = lam;
}

//=================================================================================
// main procedure

double width = 600;
double height = 600;

void SolveProblem()
{

}

void init()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, width, 0.0, height);

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
