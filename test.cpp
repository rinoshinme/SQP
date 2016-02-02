#if 0
#include <iostream>
#include "functions.h"
#include "pinv.h"
#include "pdata.h"
#include "matrix_inversion.h"
#include "utils.h"

int main()
{
	int N = 20;
	Vector x(3 * N);
	Vector mu(2);
	Vector lam(6 * N - 4);
	Vector dlax1 = dlax(x, mu, lam);

	std::getchar();
}
#endif
