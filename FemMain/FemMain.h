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
#define TagNode   1
#define TagValue  2
#define TagCheck  4
#define Static	  0
#define Model	  1
#define DynamicStatic  2
#define MSor 0
#define MLU  1
#define Explicit  0
#define Implicit  1

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
	bool * Converge;
	int SolveMethod;
	int SolveType;
	int ProblemType;
	double dT;

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

	FloatMatrix Stiff,Mass,Damp,EffictiveStiff,EffictiveMass;
	FloatArray ResultZero,LResultZero;
	FloatArray ResultFirst,LResultFirst;
	FloatArray ResultSecond,LResultSecond;
	IntArray DegreeOfFreedom;
	FloatArray ExternalForce;
	FloatArray InitialStain;
	FloatArray InteractLoad;
	FloatArray InitialDispLoad;
	FloatArray IniDisplacement;
	FloatArray InterDisplace;
	FloatArray TotalLoad;
	FloatArray EffictiveLoad;
	FloatArray Eigenvalues;

	IntArray InteractNode;
	FloatArray *InteractValue;
	FloatArray *InteractValueOld;

	LUSolve LUSolver;
	Sor Soror;

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
	void StaticSolve();
	void ModelSolve();
	void DynamicStaticSolve();
	void GloableSolve();
	bool *ConvergeCheck();
	void ComputeElementStress();
	void CountElement();
	void SendResultToNode();
	IntArray GetInteractNode(int iinter);
	FloatArray GetInteractResult(IntArray & InteractNode);
	void SetInteractResult(FloatArray * InteractResult);

	void GetSize(int & NProces);
	void GetID(int &MyID);
	FloatArray * ExchangeData();
};
