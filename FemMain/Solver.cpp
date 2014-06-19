#include "Solver.h"
#include <algorithm>

Solver::Solver()
{
}


Solver::~Solver()
{
}


int LUSolve::Decomposition(FloatMatrix &A)
{
	this->A = A;
	m = A.GetSize();
	Value = new double[m*m]();
	ipiv = new int[m]();
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Value[i*m + j] = A.at(i, j);
		}
	}
	info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, m, Value, m, ipiv);
	return info;
}
int LUSolve::Solver(FloatArray &B, FloatArray &X)
{
	this->B = B;
	this->X = X;
	if (Right == NULL)
	{
		Right = new double[m];
	}
	for (int i = 0; i < m; i++)
	{
		Right[i] = B.at(i);
	}
	trana = 'N', tranb = 'N';
	info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, trana, m, 1, Value, m, ipiv, Right, m);
	for (int i = 0; i < m; i++)
	{
		X.at(i) = Right[i];
	}
	return info;
}

bool LUSolve::Check(FloatArray &B, FloatArray &X)
{
	X.Print();
	X =A.Mult(X);
	X.Print();
	n = m; k = 1;
	lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
	if (Left == NULL)
	{
		Left = new double[m]();
	}
	for (int i = 0; i < m; i++)
	{
		Left[i] = 0;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Value[i*m + j] = A.at(i, j);
		}
	}
	A.Print();
	dgemm(&trana, &trana, &m, &n, &k, &alpha, Value, &lda, Right, &ldb, &beta, Left, &ldc);
	Error = 0;
	for (int i = 0; i < m; i++)
	{
		Error += pow((Left[i] - B.at(i)), 2);
	}
	Error = Error / m;
	if (Error < 1e-11)
	{
		return true;
	}
	return false;
}