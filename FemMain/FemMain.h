#pragma once
#pragma comment(linker,"/NODEFAULTLIB:libcmt.lib")
#ifdef WIN32
#pragma comment(lib,"GidLib\\x32\\zlib.lib")
#pragma comment(lib,"GidLib\\x32\\hdf5.lib")
#pragma comment(lib,"GidLib\\x32\\gidpost.lib")
#pragma comment(lib,"C:\\Program Files\\MPICH2\\lib\\cxx.lib")
#pragma comment(lib,"C:\\Program Files\\MPICH2\\lib\\mpi.lib")
#endif
#ifdef x64
#pragma comment(lib,"GidLib\\x64\\zlib.lib")
#pragma comment(lib,"GidLib\\x64\\hdf5.lib")
#pragma comment(lib,"GidLib\\x64\\gidpost.lib")
#pragma comment(lib,"C:\\Program Files\\MPICH2\\lib\\cxx.lib")
#pragma comment(lib,"C:\\Program Files\\MPICH2\\lib\\mpi.lib")
#endif
#include "Array.h"
#include "Model.h"
#include "Load.h"
#include "Boundary.h"
#include "Material.h"
#include "gidpost.h"
#include "Solver.h"
#include "mpi.h"

#define MAXCHAR 200

const string GroupName[4] = { " ", "Linear", " ", "Quadrilateral" };

class FemMain
{
private:
	

	string prob, text,workdir;
	string GidResultFile;
	char str[MAXCHAR];
	int nNode;
	int nElem;
	int nDim;
	int nGroup;
	int nMat;
	int nStep;
	int nDof;
	int MaxIter;
	int TotalDOF;
	int nTotalLoad;
	int nFace;
	int nConcentrate;
	int nVolumn;
	int nTotalBoundary;
	int nFix;
	int nDisp;
	int nInter;
	double Error, Tolerance;
	bool Converge;
	
	ofstream chk, res;
	ifstream glb, cor, ele, grp, loa, mat, pre;
	GiD_FILE FMesh, FRes;
	time_t tNow;
	struct tm tmLocal;
	int Id, nProces;
	
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
	FloatArray IniDisplacement;
	FloatArray InterDisplace;
	FloatArray TotalLoad;

	IntArray InteractNode;
	FloatArray InteractValue;
	FloatArray InteractValueOld;

	LUSolve *LUSolver;

public:
	void ShowTime();
	string & WorkDir();
	int GetStep();
	int GetMaxIter();
	void ReadFiles();
	void OpenGidFile();
    void GIDOutMesh();
	void GIDOutResult(int istep);
	void CloseGidFile();
	void ComputeDOF();
	void InitSolve();
	void ComputeElementStiff();
	void AssembleStiff();
	void ApplyLoad();
	void ApplyInitial();
	void Solve();
	bool ConvergeCheck();
	void ComputeElementStress();
	void CountElement();
	void SendResultToNode();
	IntArray GetInteractNode();
	FloatArray GetInteractResult(IntArray & InteractNode);
	void SetInteractResult(FloatArray &InteractResult);

	void GetSize(int & NProces);
	void GetID(int &MyID);
	FloatArray ExchangeData();
};
