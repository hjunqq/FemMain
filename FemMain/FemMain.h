#pragma once
#pragma comment(lib,"zlib.lib")
#pragma comment(lib,"hdf5.lib")
#pragma comment(lib,"gidpost.lib")
#pragma comment(linker,"/NODEFAULTLIB:libcmt.lib")
#include "Array.h"
#include "Model.h"
#include "Load.h"
#include "Boundary.h"
#include "Material.h"
#include "gidpost.h"

#define MAXCHAR 200

const string GroupName[4] = { " ", "Linear", " ", "Quadrilateral" };

class FemMain
{
private:
	

	string prob, text,workdir;
	string GidResFile;
	char str[MAXCHAR];
	int nNode;
	int nElem;
	int nDim;
	int nGroup;
	int nMat;
	int nStep;
	int nDof;
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

	FloatMatrix *Stiff;
	FloatArray *ResultZero;
	IntArray *DegreeOfFreedom;
	FloatArray *ExternalForce;
	FloatArray *InitialStain;
	FloatArray *InteractLoad;
	FloatArray *InitialDispLoad;
	FloatArray *TotalLoad;
public:
	void ShowTime();
	string & WorkDir();
	void ReadFiles();
    void GIDOutMesh();
	void ComputeDOF();
	void InitSolve();
	void ComputeElementStiff();
};
