#pragma once
#include "Array.h"
#include "Model.h"
#include "Load.h"
#include "Boundary.h"
#include "Material.h"
class FemMain
{
private:
	string prob, text;
	string GidResFile;
	int nNode;
	int nElem;
	int nDim;
	int nGroup;
	int nMat;
	int nStep;
	int TotalDOF;
	
	ofstream chk, res;
	time_t tNow;
	struct tm tmLocal;
	int Id, nPorcs;
	
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
};
void FemMain::ShowTime()
{
	time(&tNow);
	localtime_s(&tmLocal, &tNow);
	cout << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon + 1 << "-" <<
		tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0') << setw(2) << tmLocal.tm_min << ":"
		<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
	chk << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon << "-" <<
		tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0') << setw(2) << tmLocal.tm_min << ":"
		<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
}
