#pragma once
#include "Array.h"
#include "mkl.h"
class Solver
{
public:
	Solver();
	virtual ~Solver();
};
class LUSolve
{
private:
	FloatMatrix *A;
	FloatArray *X;
	FloatArray *B;
	int *ipiv;
public:
	int Decomposition(FloatMatrix *A);
	int Solver(FloatMatrix *A, FloatArray *B, FloatArray *X);
};

class ElasticSolver :
	public Solver
{
private:
	double ErrorLimit;
	int MaxIterate;

};

