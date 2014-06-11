#pragma once
#pragma comment(lib,"zlib.lib")
#pragma comment(lib,"hdf5.lib")
#pragma comment(lib,"gidpost.lib")
#include "Array.h"
#include "Model.h"
#include "Load.h"
#include "Boundary.h"
#include "Material.h"
#include "gidpost.h"

#define MAXCHAR 200
class FemMain
{
private:
	string GroupName[4] = {"", "Linear","", "Quadrilateral" };

	string prob, text,workdir;
	string GidResFile;
	char str[MAXCHAR];
	int nNode;
	int nElem;
	int nDim;
	int nGroup;
	int nMat;
	int nStep;
	int TotalDOF;
	int nTotalLoad;
	int nFace;
	int nConcentrate;
	int nVolumn;
	int nTotalBoundary;
	int nFix;
	int nDisp;
	int nInter;
	
	ofstream chk, res;
	ifstream glb, cor, ele, grp, loa, mat, pre;
	time_t tNow;
	struct tm tmLocal;
	int Id, nPorcs;
	
	Group *Groups;
	Material *Mats;
	Node *Nodes;
	Face *Faces;
	Concentrate *Cons;
	Volumn *Vols;
	Fixed *Fix;
	Displace *Disp;
	Interact *Inters;

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
    void GIDOutMesh();
};
