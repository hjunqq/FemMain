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
class SOR
{

};
class ElasticSolver :
	public Solver
{
private:
	double ErrorLimit;
	int MaxIterate;

};

