#include "utils.h"
#include <stdexcept>
#include <iostream>

#if 0
double dist(double x1, double y1, double x2, double y2)
{
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

double dist2(double x1, double y1, double x2, double y2)
{
	return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
}
#endif

void swap(double* a, double* b)
{
	double c;
	c = *a;
	*a = *b;
	*b = c;
}

int inv(double* p, int n)
{
	int *is, *js;
	int i, j, k;

#if 0
	// print matrix;
	for (i = 0; i < n; ++i)
	{
		putchar('\n');
		for (j = 0; j < n; ++j)
			printf("%f   ", *(p + i * n + j));
	}
	puts("\n\n\n\n");
#endif

	double temp, fmax;
	is = (int*)malloc(n * sizeof(int));
	js = (int*)malloc(n * sizeof(int));

	for (k = 0; k < n; ++k)
	{
		fmax = 0.0;
		for (i = k; i < n; ++i)
		{
			for (j = k; j < n; ++j)
			{
				temp = fabs(*(p + i*n + j));
				if (temp > fmax)
				{
					fmax = temp;
					is[k] = i;
					js[k] = j;
				}
			}
		}
		if ((fmax + 1.0) == 1.0)
		{
			free(is);
			free(js);
			return 0;
		}

		if ((i = is[k]) != k)
			for (j = 0; j < n; ++j)
				swap(p+k*n+j, p+i*n+j);
		if ((j = js[k]) != k)
			for (i = 0; i < n; ++i)
				swap(p+i*n+k, p+i*n+j);
		p[k*n+k] = 1.0/p[k*n+k];

		for (j = 0; j < n; ++j)
			if (j != k)
				p[k*n+j] *= p[k*n+k];
		for (i = 0; i < n; ++i)
			if (i != k)
				for (j = 0; j < n; ++j)
					if (j != k)
						p[i*n+j] = p[i*n+j] - p[i*n+k]*p[k*n+j];
		for (i = 0; i < n; ++i)
			if (i != k)
				p[i*n+k] *= -p[k*n+k];
	}
	for (k = n-1; k >= 0; --k)
	{
		if ((j = js[k]) != k)
			for (i = 0; i < n; ++i)
				swap(p+j*n+i, p+k*n+i);
		if ((i=is[k]) != k)
			for (j = 0; j < n; ++j)
				swap(p+j*n+i, p+j*n+k);
	}
	free(is);
	free(js);
	return 1;
}


// ===============================================================================
// Vector operators

Vector& Vector::operator+=(const Vector& rhs)
{
	for (int i = 0; i < GetSize(); ++i)
		x[i] += rhs.x[i];
	return *this;
}

Vector& Vector::operator-=(const Vector& rhs)
{
	for (int i = 0; i < GetSize(); ++i)
		x[i] -= rhs.x[i];
	return *this;
}

Vector operator+(const Vector& lhs, const Vector& rhs)
{
	Vector result(lhs);
	return result += rhs;
}

Vector operator-(const Vector& lhs, const Vector& rhs)
{
	Vector result(lhs);
	return result-=rhs;
}

Vector& Vector::operator*=(double a)
{
	for (int i = 0; i < GetSize(); ++i)
		x[i] *= a;
	return *this;
}

Vector& Vector::operator/=(double a)
{
	for (int i = 0; i < GetSize(); ++i)
		x[i] /= a;
	return *this;
}

Vector operator*(const Vector& vector, double a)
{
	Vector result(vector);
	return result *= a;
}

Vector operator/(const Vector& vector, double a)
{
	Vector result(vector);
	return result /= a;
}

Vector operator*(double a, const Vector& vector)
{
	Vector result(vector);
	return result *= a;
}

// =======================================================================================
// Matrix operators

Matrix& Matrix::operator+=(const Matrix& rhs)
{
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			SetValue(i, j, GetValue(i, j) + rhs.GetValue(i, j));
		}
	}
	return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs)
{
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			SetValue(i, j, GetValue(i, j) - rhs.GetValue(i, j));
		}
	}
	return *this;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
	Matrix result(lhs);
	return result += rhs;
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
	Matrix result(lhs);
	return result -= rhs;
}

Matrix& Matrix::operator*=(double a)
{
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			SetValue(i, j, GetValue(i, j) * a);
		}
	}
	return *this;
}

Matrix& Matrix::operator/=(double a)
{
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			SetValue(i, j, GetValue(i, j) / a);
		}
	}
	return *this;
}

Matrix operator*(const Matrix& matrix, double a)
{
	Matrix result(matrix);
	return result *= a;
}

Matrix operator*(double a, const Matrix& matrix)
{
	return matrix * a;
}

Matrix operator/(const Matrix& matrix, double a)
{
	Matrix result(matrix);
	return result /= a;
}

// =====================================================================================
// Matrix and Vector manipulations

Vector VxM(const Vector& v, const Matrix& m)
{
	int row = m.GetNumRows();
	int col = m.GetNumCols();
	Vector result(col);
	for (int i = 0; i < col; ++i)
	{
		double sum = 0.0;
		for (int k = 0; k < row; ++k)
			sum += v.GetValue(k) * m.GetValue(k, i);
		result.SetValue(i, sum);
	}
	return result;
}

Vector MxV(const Matrix& m, const Vector& v)
{
	int row = m.GetNumRows();
	int col = m.GetNumCols();

	Vector result(row);
	for (int i = 0; i < row; ++i)
	{
		double sum = 0.0;
		for (int k = 0; k < col; ++k)
			sum += m.GetValue(i, k) * v.GetValue(k);
		result.SetValue(i, sum);
	}
	return result;
}

void MxM(const Matrix& m1, const Matrix& m2, Matrix* result)
{
	int row1 = m1.GetNumRows();
	int col1 = m1.GetNumCols();
	int row2 = m2.GetNumRows();
	int col2 = m2.GetNumCols();

	if (col1 != row2)
		return;

	result->Resize(row1, col2);
	for (int i = 0; i < row1; ++i)
	{
		for (int j = 0; j < col2; ++j)
		{
			double sum = 0.0;
			for (int k = 0; k < col1; ++k)
				sum += m1.GetValue(i, k) * m2.GetValue(k, j);
			result->SetValue(i, j, sum);
		}
	}
}

Matrix VStack(const Matrix& m1, const Matrix& m2)
{
	int row1 = m1.GetNumRows();
	int row2 = m2.GetNumRows();
	int col1 = m1.GetNumCols();
	int col2 = m2.GetNumCols();
	if (col1 != col2)
		throw std::runtime_error("Col size not equal");
	Matrix result(row1+row2, col1);
	for (int i = 0; i < row1; ++i)
		for (int j = 0; j < col1; ++j)
			result.SetValue(i, j, m1.GetValue(i, j));
	for (int i = 0; i < row2; ++i)
		for (int j = 0; j < col2; ++j)
			result.SetValue(i+row1, j, m2.GetValue(i, j));
	return result;
}

double Dot(const Vector& v1, const Vector& v2)
{
	double sum = 0.0;
	int n = v1.GetSize();
	for (int i = 0; i < n; ++i)
		sum += v1.GetValue(i) * v2.GetValue(i);
	return sum;
}

// outer product -> matrix
Matrix Out(const Vector& v1, const Vector& v2)
{
	int row = v1.GetSize();
	int col = v2.GetSize();

	Matrix result(row, col);
	for (int i = 0; i < row; ++i)
		for (int j = 0; j < col; ++j)
			result.SetValue(i, j, v1.GetValue(i) * v2.GetValue(j));
	return result;
}

double vBv(const Vector& v, const Matrix& B)
{
	Vector temp = VxM(v, B);
	return Dot(temp, v);
}

void DisplayVector(const Vector& v)
{
	for (int i = 0; i < v.GetSize(); ++i)
		std::cout << v.GetValue(i) << std::endl;
	std::cout << std::endl;
}

void DisplayMatrix(const Matrix& m)
{
	int rows = m.GetNumRows();
	int cols = m.GetNumCols();
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			std::cout << m.GetValue(i, j) << "   ";
		}
		std::cout << std::endl;
	}
}

Vector LeastSquare(const Matrix& matrix, const Vector& vector)
{
	// 1. construct Eigen matrix and Eigen vector
	// 2. solve the least square problem 
	// 3. transform back to Vector object
	int row = matrix.GetNumRows();
	int col = matrix.GetNumCols();
	int len = vector.GetLength();

	Eigen::MatrixXd mat(col, row);
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
			mat(j, i) = matrix.GetValue(i, j);
	}

	Eigen::VectorXd vec(len);
	for (int i = 0; i < len; ++i)
		vec(i) = vector.GetValue(i);

	Eigen::VectorXd res = mat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vec);

	Vector result(row);
	for (int i = 0; i < row; ++i)
		result.SetValue(i, res(i));
	return result;
}
