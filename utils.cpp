#include "utils.h"
#include <stdexcept>

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
		for (int j = 0; j < col1; ++j)
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
