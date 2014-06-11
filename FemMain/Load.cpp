#include "Load.h"


Load::Load()
{
	Index = 0;
}


Load::~Load()
{
}


Face::Face()
{
	nElem = 0;
	StratV = 0.0;
	EndV = 0.0;
	Dir = 0;
	StartC = 0.0;
	EndC = 0.0;
}


Face::~Face()
{
}


void Face::Init(int Index, int nEle, int nNode, double StartC, double StartV, double EndC, double EndV, int Dir)
{
	this->Index = Index;
	this->nElem = nEle;
	this->nNode = nNode;
	this->StartC = StartC;
	this->StratV = StartV;
	this->EndC = EndC;
	this->EndV = EndV;
	this->Dir = Dir;
	if (nNode == 2)
	{
		Lines = new Line[nElem]();
	}
}
void Face::Set(Line *TLines)
{
	for (int i = 0; i < nElem; i++)
	{
		Lines[i] =TLines[i];
	}
}


Volumn::Volumn()
{
	Index = 0;
	group = 0;
	Dir = 0;
	Value = 0.0;
}


Volumn::~Volumn()
{
}


void Volumn::Init(int Index, int nGroup,IntArray *group, double Value, int Dir)
{
	this->Index = Index;
	this->Value = Value;
	this->Dir = Dir;
	this->nGroup = nGroup;
	this->group = new IntArray(*group);
}


Concentrate::Concentrate()
{
	nNode = 0;
	Dir = 0;
}


Concentrate::~Concentrate()
{
}


void Concentrate::Init(int Index,int nNode, IntArray * Nodes, FloatArray * Values, int Dir)
{
	this->Index = Index;
	this->nNode = nNode;
	this->Dir = Dir;
	this->Nodes = new IntArray(*Nodes);
	this->Values = new FloatArray(*Values);
}
