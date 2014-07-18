#pragma once
#include "Array.h"
#include "mkl.h"
#include <algorithm>
class Solver
{
public:
	Solver();
	virtual ~Solver();
};
class LUSolve
{
private:
	int m,n,info,k,lda,ldb,ldc;
	char trana, tranb;
	double *Value, *Right,*Left;
	FloatMatrix A;
	FloatArray X;
	FloatArray B;
	int *ipiv;
	double alpha = 1.0, beta = 1.0;
	double Error;
public:
	int Decomposition(FloatMatrix &A);
	int Inverse();
	void Mult(FloatArray &B);
	int Solver(FloatArray &B, FloatArray &X);
	bool Check(FloatArray &B, FloatArray &X);
};
class Sor
{
private:
	int m;
	double w,SpectralRadius;
	double axnn, axn;
	double *Value,*X1,*X2;
	FloatMatrix A;
	FloatArray V1,V2;
	double Error,ErrorAverage;
public:
	void Init(FloatMatrix &A);
	void MaxSpectralRadius();
	void Compute(FloatArray &B);
	void Solve(FloatArray &B,FloatArray &X);
};
class Newmark:
	public Solver
{
private:
	double dT;
	double Alpha;
	double Beta;
	double a0, a1, a2, a3, a4, a5, a6, a7;
	FloatMatrix A;
	FloatArray X;
	FloatArray B;
	int *ipiv;
	double Error;
public:
	int IntSolver(double Alpha, double Beta);
	int SetStiffMatix(FloatMatrix &A);

};
class ElasticSolver :
	public Solver
{
private:
	double ErrorLimit;
	int MaxIterate;

};

