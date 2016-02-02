#include "pinv.h"

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
