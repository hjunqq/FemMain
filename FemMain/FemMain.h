#pragma once
#pragma comment(linker,"/NODEFAULTLIB:libcmt.lib")
#ifdef WIN32
#pragma comment(lib,"GidLib\\x32\\zlib.lib")
#pragma comment(lib,"GidLib\\x32\\hdf5.lib")
#pragma comment(lib,"GidLib\\x32\\gidpost.lib")
#endif
#ifdef x64
#pragma comment(lib,"GidLib\x64\zlib.lib")
#pragma comment(lib,"GidLib\x64\hdf5.lib")
#pragma comment(lib,"GidLib\x64\gidpost.lib")
#endif
#include "Array.h"
#include "Model.h"
#include "Load.h"
#include "Boundary.h"
#include "Material.h"
#include "gidpost.h"
#include "Solver.h"

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

	FloatMatrix Stiff;
	FloatArray ResultZero;
	IntArray DegreeOfFreedom;
	FloatArray ExternalForce;
	FloatArray InitialStain;
	FloatArray InteractLoad;
	FloatArray InitialDispLoad;
	FloatArray TotalLoad;

	LUSolve *LUSolver;

public:
	void ShowTime();
	string & WorkDir();
	int GetStep();
	void ReadFiles();
    void GIDOutMesh();
	void GIDOutResult();
	void ComputeDOF();
	void InitSolve();
	void ComputeElementStiff();
	void AssembleStiff();
	void ApplyLoad();
	void ApplyInitial();
	void Solve();
	bool ConvergeCheck();
	void ComputeElementStress();
};
