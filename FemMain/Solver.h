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
	FloatArray V1,V2,U1,U2;
	double Error,ErrorAverage;
public:
	void Init(FloatMatrix &A);
	void MaxSpectralRadius();
	void Compute(FloatArray &B,FloatArray &X);
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
	FloatArray m1, m2, m3;
	FloatArray c1, c2, c3;
	FloatArray MassArray, DampArray;
	FloatArray Acc1, Acc2, Acc3;
	FloatMatrix StiffEffictive;
	FloatArray LoadEfficitive;
	FloatArray X;
	FloatArray B;
	int *ipiv;
	double Error;
public:
	void IntSolver(double dT);
	FloatMatrix & EffictiveStiff(FloatMatrix &Stiff, FloatMatrix &Mass, FloatMatrix &Damp);
	FloatArray & EffictiveLoad(FloatArray &Load, FloatArray &ResultZero, FloatArray &ResultFirst, FloatArray &ResultSecond,
		const FloatMatrix Mass, const FloatMatrix Damp);
	void SolvePorcess(const FloatArray ResultZero,const FloatArray LResultZero,
		FloatArray &ResultFirst,const FloatArray LResultFirst,
		FloatArray &ResultSecond,const FloatArray LResultSecond);
};
class ElasticSolver :
	public Solver
{
private:
	double ErrorLimit;
	int MaxIterate;

};
class CentralDifference:
	public Solver
{
private:
	double dT;
	double c0, c1, c2, c3;
	FloatMatrix MassEffictive;
	FloatArray LoadEffictive;
	FloatArray LastDisplace,InitAcc;
	FloatArray Acc, Velocity;
public:
	void Init(double dT, int TotDOF);
	FloatArray & Init(FloatArray &ResultSecond);
	FloatMatrix & EffictiveMass(FloatMatrix & Mass, FloatMatrix & Damp);
	FloatArray &EffictiveLoad(FloatArray &Load, FloatMatrix & Stiff, FloatMatrix &Mass,
		FloatMatrix & Damp, FloatArray ResultZero, FloatArray LResultZero);
	void SolvePorcess(const FloatArray ResultZero, const FloatArray LResultZero,
		FloatArray &ResultFirst, FloatArray &ResultSecond);
};