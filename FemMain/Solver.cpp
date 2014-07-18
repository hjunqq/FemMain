#include "Solver.h"


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
int LUSolve::Inverse()
{
	info = LAPACKE_dgetri(LAPACK_COL_MAJOR, m, Value, m, ipiv);
	return info;
}
void LUSolve::Mult(FloatArray &B)
{
	k = m; n = 1;
	lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
	if (Left == NULL)
	{
		Left = new double[m]();
	}
	for (int i = 0; i < m; i++)
	{
		Left[i] = 0;
	}
	if (Right == NULL)
	{
		Right = new double[m];
	}
	for (int i = 0; i < m; i++)
	{
		Right[i] = B.at(i);
	}
	trana = 'N', tranb = 'N';
	dgemm(&trana, &trana, &m, &n, &k, &alpha, Value, &lda, Right, &ldb, &beta, Left, &ldc);
	for (int i = 0; i < m; i++)
	{
		B.at(i) = Left[i];
	}
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
	k = m; n = 1;
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

void Sor::Init(FloatMatrix &A)
{
	this->A = A;
	m = A.GetSize();
	Value = new double [m*m]();
	X1 = new double[m]();
	X2 = new double[m]();
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Value[i*m + j] = A.at(i, j);
		}
	}
	for (int i = 0; i < m; i++)
	{
		X1[i] = 1;
		X2[i] = 0;
	}
	MaxSpectralRadius();
	w = 2 / (1 + sqrt(1 - SpectralRadius*SpectralRadius));
}
void Sor::MaxSpectralRadius()
{
	V1.SetSize(m);
	V2.SetSize(m);
	V1.Set(1);
	FloatMatrix D, LU,J;
	D.SetSize(m, m);
	LU = A;
	for (int i = 0; i < m; i++)
	{
		D.at(i, i) = 1/A.at(i, i);
		LU.at(i, i) = 0;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			LU.at(i, j) = -LU.at(i, j);
		}
	}
	J = D.Mult(LU);
	double OldSpectralRadius = 0, Error = 0;
	do
	{
		V2 = J.Mult(V1);
		SpectralRadius = V2.at(0) / V1.at(0);
		Error = abs(SpectralRadius - OldSpectralRadius);
		OldSpectralRadius = SpectralRadius;
	} while (Error > 10e-5);
}
void Sor::Compute(FloatArray &B)
{
	for (int i = 0; i < m; i++)
	{
		axnn = axn = 0;
		for (int j = 0; j < i - 1; j++)
		{
			axnn += Value[i*m + j] * X2[j];
		}
		for (int j = 0; j < m; j++)
		{
			axn += Value[i*m + j] * X1[j];
		}
		X2[i] = X1[i] + w*(B.at(i)-axnn-axn) / Value[i*m + i];
	}
}

void Sor::Solve(FloatArray &B,FloatArray &X)
{
	//ErrorAverage = 0;
	//for (int i = 0; i < m; i++)
	//{
	//	ErrorAverage += B.at(i)*B.at(i);
	//}
	//ErrorAverage = ErrorAverage / m;
	//if (0 == ErrorAverage )
	//{
	//	ErrorAverage = 1;
	//}
	int iiter = 0;
	do
	{
		iiter++;
		Compute(B);
		Error = 0;
		for (int i = 0; i < m; i++)
		{
			Error += abs(X2[i] - X1[i]);
		}
		for (int i = 0; i < m; i++)
		{
			X1[i] = X2[i];
		}
		Error = Error / m;
		cout << "iiter=       " << iiter << "    Error=      " << Error << endl;
	} while (Error > 10e-5);
	for (int i = 0; i < m; i++)
	{
		X.at(i) = X2[i];
	}
}