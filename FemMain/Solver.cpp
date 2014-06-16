#include "Solver.h"


Solver::Solver()
{
}


Solver::~Solver()
{
}


int LUSolve::Decomposition(FloatMatrix *A)
{
	int m;
	m = A->GetSize();
	double *Value;
	Value = new double[m*m]();
	ipiv = new int[m]();
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Value[i*m + j] = A->at(i, j);
		}
	}
	int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, m, Value, m, ipiv);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			A->at(i, j) = Value[i*m + j];
		}
	}
	A->Print();
	return info;
}
int LUSolve::Solver(FloatMatrix *A, FloatArray *B, FloatArray *X)
{
	int m;
	char trana;
	double *Value,*Right;
	m = A->GetSize();
	Value = new double[m*m]();
	Right = new double[m]();
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Value[i*m + j] = A->at(i, j);
		}
	}
	for (int i = 0; i < m; i++)
	{
		Right[i] = B->at(i);
	}

	trana = 'N';
	
	int info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, trana, m, 1, Value, m, ipiv, Right, m);
	for (int i = 0; i < m; i++)
	{
		X->at(i) = Right[i];
	}

	return info;
}