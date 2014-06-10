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


void Face::Init(int Index, double StartC, double StartV, double EndC, double EndV, int Dir)
{
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


void Volumn::Init(int Index, int group, double Value, int Dir)
{
}


Concentrate::Concentrate()
{
	nNode = 0;
	Dir = 0;
}


Concentrate::~Concentrate()
{
}


void Concentrate::Init(int Index, IntArray * Nodes, FloatArray * Values, int Dir)
{
}
