#pragma once
#include "Array.h"
#include "Model.h"
#include "Load.h"
#include "Boundary.h"
#include "Material.h"
class FemMain
{
private:
	string prob, text,workdir;
	string GidResFile;
	int nNode;
	int nElem;
	int nDim;
	int nGroup;
	int nMat;
	int nStep;
	int TotalDOF;
	
	ofstream chk, res;
	ifstream glb, cor, ele, grp, loa, mat, pre;
	time_t tNow;
	struct tm tmLocal;
	int Id, nPorcs;
	
	Group *Groups;
	Quadr *Quadrs;
	Material *Mats;
	Node *Nodes;

	FloatMatrix *A;
	FloatArray *Result;
	IntArray DegreeOfFreedom;
	FloatArray *LoadMatrix;
	FloatArray *InitialLoad;
	FloatArray *InteractLoad;
	FloatArray *DispLoad;
	FloatArray *TotalLoad;
public:
	void ShowTime();
	string & WordDir();
	void ReadFiles();
};
