#ifndef _SQP_UTILS_H
#define _SQP_UTILS_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <memory>
#include "matrix_inversion.h"

class Vector
{
private:
	double* x;
	int length;

public:
	/* constructors and destructor */
	Vector()
		:length(0), x(NULL)
	{}

	Vector(double* data, int n)
		:length(n)
	{
		x = new double[length];
		memcpy(x, data, sizeof(double) * length);
	}

	Vector(int len)
		:length(len)
	{
		x = new double[length];
		memset(x, 0, sizeof(double)*length);
	}

	Vector(const Vector& copy)
	{
		length = copy.length;
		x = new double[length];
		memcpy(x, copy.x, sizeof(double) * length);
	}

	Vector& operator=(const Vector& rhs)
	{
		if (this == &rhs)
			return *this;
		length = rhs.length;
		if (x)
			delete[] x;
		x = new double[length];
		memcpy(x, rhs.x, sizeof(double) * length);
		return *this;
	}

	double& operator[](int index)
	{
		return x[index];
	}

	~Vector()
	{
		if (x)
			delete[] x;
	}

	// methods
	int GetLength() const
	{
		return length;
	}

	int GetSize() const
	{
		return length;
	}

	double* GetData() const
	{
		return x;
	}

	double GetValue(int index) const
	{
		return x[index];
	}

	void SetValue(int index, double value)
	{
		x[index] = value;
	}

	double Norm() const
	{
		double sum = 0.0;
		for (int i = 0; i < GetSize(); ++i)
			sum += x[i] * x[i];
		return sqrt(sum);
	}

	// sum over all elements
	double Norm1() const
	{
		double sum = 0.0;
		for (int i = 0; i < GetSize(); ++i)
			sum += x[i];
		return sum;
	}

	// find abs maximum value
	double NormInf() const
	{
		double max = 0.0;
		for (int i = 0; i < GetSize(); ++i)
		{
			if (fabs(x[i]) > max)
				max = fabs(x[i]);
		}
		return max;
	}

	Vector& operator+=(const Vector& rhs);
	Vector& operator-=(const Vector& rhs);
	friend Vector operator+(const Vector& lhs, const Vector& rhs);
	friend Vector operator-(const Vector& lhs, const Vector& rhs);

	Vector& operator*=(double a);
	Vector& operator/=(double a);
	friend Vector operator*(const Vector& vector, double a);
	friend Vector operator/(const Vector& vector, double a);
	friend Vector operator*(double a, const Vector& vector);
};

class Matrix
{
private:
	double* data;
	int rows;
	int cols;

public:
	/* constructors and destructor */
	Matrix()
		:data(NULL), rows(0), cols(0)
	{}

	Matrix(int row, int col, double* d=NULL)
		:rows(row), cols(col)
	{
		data = new double[rows * cols];
		if (!d)
			memset(data, 0, sizeof(double) * rows * cols);
		else
			memcpy(data, d, sizeof(double) * rows * cols);
	}

	Matrix& operator=(const Matrix& rhs)
	{
		if (this == &rhs)
			return *this;
		rows = rhs.rows;
		cols = rhs.cols;
		if (data)
			delete[] data;
		data = new double[rows * cols];
		memcpy(data, rhs.data, sizeof(double) * rows * cols);
		return *this;
	}

	~Matrix()
	{
		if (data)
			delete[] data;
	}

	// methods
	int GetNumRows() const
	{
		return rows;
	}

	int GetNumCols() const
	{
		return cols;
	}

	double GetValue(int row, int col) const
	{
		return data[row * cols + col];
	}

	void SetValue(int row, int col, double value)
	{
		data[row * cols + col] = value;
	}

	Matrix Transpose() const
	{
		int newRow = cols;
		int newCol = rows;
		Matrix result(newRow, newCol);
		for (int i = 0; i < newRow; ++i)
			for (int j = 0; j < newCol; ++j)
				result.SetValue(i, j, GetValue(j, i));
		return result;
	}

	void Resize(int r, int c)
	{
		rows = r;
		cols = c;
		if (data)
			delete[] data;
		data = new double[rows * cols];
		memset(data, 0, sizeof(double) * rows * cols);
	}

	void Reset()
	{
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				SetValue(i, j, 0.0);
	}

	// check rows == cols
	Matrix Inverse() const
	{
		double* m = new double[rows * cols];
		memcpy(m, data, sizeof(double) * rows * cols);
		inv(m, rows);
		return Matrix(rows, rows, m);
	}

	Matrix& operator+=(const Matrix& rhs);
	Matrix& operator-=(const Matrix& rhs);
	friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
	friend Matrix operator-(const Matrix& lhs, const Matrix& rhs);

	Matrix& operator*=(double a);
	Matrix& operator/=(double a);
	friend Matrix operator*(const Matrix& matrix, double a);
	friend Matrix operator*(double a, const Matrix& matrix);
	friend Matrix operator/(const Matrix& matrix, double a);
	
};

// Matrix and Vector manipulations
Vector VxM(const Vector& v, const Matrix& m);
Vector MxV(const Matrix& m, const Vector& v);
void MxM(const Matrix& m1, const Matrix& m2, Matrix* result);

Matrix VStack(const Matrix& m1, const Matrix& m2);
double Dot(const Vector& v1, const Vector& v2); // inner product -> double
Matrix Out(const Vector& v1, const Vector& v2); // outer product -> matrix
double vBv(const Vector& v, const Matrix& B);

#ifdef __cplusplus
extern "C" {
#endif

double dist(double x1, double y1, double x2, double y2);

double dist2(double x1, double y1, double x2, double y2);

#ifdef __cplusplus
}
#endif
#endif

