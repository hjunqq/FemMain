#include "Solver.h"


Solver::Solver()
{
}


Solver::~Solver()
{
}

LUSolve::LUSolve()
{
	Value = NULL;
	ipiv = NULL;
	Right = NULL;
	Left = NULL;
	Error = 0;
	info = 0;
	m = n = k = 0;
	lda = ldb = ldc = 0;
	trana =tranb ='N';
}
LUSolve::~LUSolve()
{
	if (Value != NULL)
	{
		delete[] Value;
		Value = NULL;
	}
	if (ipiv != NULL)
	{
		delete[] ipiv;
		ipiv = NULL;
	}
}
int LUSolve::Decomposition(FloatMatrix &A)
{
	this->A = A;
	m = A.GetSize();
	if (Value != NULL)
	{
		delete[] Value;
		Value = NULL;
	}
	Value = new double[m*m]();
	if (ipiv != NULL)
	{
		delete[] ipiv;
		ipiv = NULL;
	}
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
Sor::Sor()
{
	Value = NULL;
	X1 = NULL;
	X2=NULL;
	axnn = 0;
	axn = 0;
	Error = 0;
	m = 0;
	w = 0;
	SpectralRadius = 0;
	ErrorAverage = 0;
}
Sor::~Sor()
{
	if (Value != NULL)
	{
		delete[] Value;
		Value = NULL;
	}
	if (X1 != NULL)
	{
		delete[] X1;
		X1 = NULL;
	}
	if (X2 != NULL)
	{
		delete[]X2;
		X2 = NULL;
	}
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
	w = 1.2;
}
void Sor::MaxSpectralRadius()
{
	V1.SetSize(m);
	V2.SetSize(m);
	U1.SetSize(m);
	U2.SetSize(m);
	V1.Set(1);
	U1.Set(1);
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
	int iiter = 0;
	do
	{
		iiter++;
		V1 = J.Mult(U1);
		SpectralRadius = V1.MaxValue();
		U1 = V1.Times(1 / SpectralRadius);
		Error = abs(SpectralRadius - OldSpectralRadius);
		OldSpectralRadius = SpectralRadius;		
	} while (iiter<20 && Error > 10e-5);
}
void Sor::Compute(FloatArray &B, FloatArray &X)
{
	for (int i = 0; i < m; i++)
	{
		axnn = axn = 0;
		for (int j = 0; j < i ; j++)
		{
			axnn += Value[i*m + j] * X2[j];
		}
		for (int j = i+1 ; j < m; j++)
		{
			axn += Value[i*m + j] * X1[j];
		}
		X2[i] = (1-w)*X1[i] + w*(B.at(i)-axnn-axn) / Value[i*m + i];
		_ASSERT(Value[i*m + i] != 0);
	}
	for (int i = 0; i < m; i++)
	{
		X.at(i) = X2[i];
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
	FloatArray T;
	T.SetSize(m);
	int iiter = 0;
	do
	{
		iiter++;
		Compute(B,T);
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
	} while ( Error > 10e-11);
	for (int i = 0; i < m; i++)
	{
		X.at(i) = X2[i];
	}
}
Newmark::Newmark()
{
	ipiv = NULL;
	a0 = a1 = a2 = a3 = a4 = a5 = a6 = a7 =0;
	Alpha = 0;
	Beta = 0;
	Error = 0;
	dT = 0;
}
Newmark::~Newmark()
{
	if (ipiv != NULL)
	{
		delete[] ipiv;
		ipiv = NULL;
	}
}
void Newmark::IntSolver(double dT)
{
	Alpha = 0.25; Beta = 0.5;
	a0 = 1 / (Alpha*dT*dT);
	a1 = Beta / (Alpha*dT);
	a2 = 1 / (Alpha*dT);
	a3 = 1 / (2 * Alpha) - 1;
	a4 = Beta / Alpha - 1;
	a5 = dT / 2 * (Beta / Alpha - 2);
	a6 = (1 - Beta)*dT;
	a7 = Beta*dT;
}

FloatMatrix & Newmark::EffictiveStiff(FloatMatrix &Stiff, FloatMatrix &Mass, FloatMatrix &Damp)
{
	StiffEffictive = Stiff + Mass.Mult(a0) + Damp.Mult(a1);
	return StiffEffictive;
}

FloatArray &Newmark::EffictiveLoad(FloatArray &Load, FloatArray &ResultZero, FloatArray &ResultFirst, FloatArray &ResultSecond,
	const FloatMatrix Mass,const FloatMatrix Damp)
{
	m1 = ResultZero.Times(a0);
	m2 = ResultFirst.Times(a2);
	m3 = ResultSecond.Times(a3);
	c1 = ResultZero.Times(a1);
	c2 = ResultFirst.Times(a4);
	c3 = ResultSecond.Times(a5);
	MassArray = m1 + m2 + m3;
	DampArray = c1 + c2 + c3;
	MassArray = Mass.Mult(MassArray);
	DampArray = Damp.Mult(DampArray);
	LoadEfficitive = Load + MassArray +DampArray ;
	return LoadEfficitive;
}
void Newmark::SolvePorcess(const FloatArray ResultZero, const FloatArray LResultZero,
	FloatArray &ResultFirst, const FloatArray LResultFirst,
	FloatArray &ResultSecond, const FloatArray LResultSecond)
{
	Acc1 = ResultZero - LResultZero;
	Acc1=Acc1.Times(a0);
	Acc2 = ResultFirst.Times(a2);
	Acc3 = ResultSecond.Times(a3);
	ResultSecond = Acc1 - Acc2 - Acc3;
	ResultFirst = LResultFirst + LResultSecond.Times(a6) + ResultSecond.Times(a7);
}
CentralDifference::CentralDifference()
{
	c0 = c1 = c2 = c3 = 0;
	dT = 0;
}
CentralDifference::~CentralDifference()
{

}
void CentralDifference::Init(double dT,  int TotDOF)
{
	this->dT = dT;
	c0 = 1 / (dT*dT);
	c1 = 1 / (2 * dT);
	c2 = 2 * c0;
	c3 = 1 / c2;
	LastDisplace.SetSize(0);
}
FloatArray & CentralDifference::Init(FloatArray &ResultSecond)
{
	InitAcc = ResultSecond;
	LastDisplace = InitAcc.Times(c3);
	return LastDisplace;
}
FloatMatrix & CentralDifference::EffictiveMass(FloatMatrix & Mass, FloatMatrix & Damp)
{
	MassEffictive = Mass.Mult(c0)+Damp.Mult(c1);
	return MassEffictive;
}
FloatArray & CentralDifference::EffictiveLoad(FloatArray &Load, FloatMatrix & Stiff, FloatMatrix &Mass,
	FloatMatrix & Damp, FloatArray ResultZero, FloatArray LResultZero)
{
	FloatMatrix KM, M,C;
	FloatArray at, lastat; 
	LoadEffictive = Load;
	KM = Mass.Mult(c2);
	KM = Stiff - KM;
	at = KM.Mult(ResultZero);
	
	M = Mass.Mult(c0);
	C = Damp.Mult(c1);
	M = M - C;
	lastat = M.Mult(LResultZero);
	//cout << endl;
	//cout << endl;
	//cout << "Load        ";
	//LoadEffictive.Print();
	//cout << "at          ";
	//at.Print();
	//cout << "LastAt      ";
	//lastat.Print();
	LoadEffictive = LoadEffictive - at - lastat;
	return LoadEffictive;
}
void CentralDifference::SolvePorcess(const FloatArray ResultZero, const FloatArray LResultZero,
	FloatArray &ResultFirst,FloatArray &ResultSecond)
{
	Acc = LastDisplace - LResultZero.Times(2) - ResultZero;
	ResultSecond = Acc.Times(c0);
	Velocity = LastDisplace.Times(-1) + ResultZero;
	ResultSecond = Velocity.Times(c1);
	
}
