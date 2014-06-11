#pragma once
#include "Array.h"
#include "Model.h"
class Load
{
public:
	Load();
	virtual ~Load();
protected:
	int Index;
	Line *Lines;
	Quadr *Quadrs;
};

class Face :
	public Load
{
public:
	Face();
	virtual ~Face();
private:
	int nElem;
	int nNode;
	double StratV;
	double EndV;
	int Dir;
	double StartC;
	double EndC;
public:
	void Init(int Index,int nEle,int nNode, double StartC, double StartV, double EndC, double EndV, int Dir);
	void Set(Line *Lines);
	void Set(Quadr *Quadrs);
};

class Volumn :
	public Load
{
public:
	Volumn();
	virtual ~Volumn();
private:
	int Index;
	int nGroup;
	IntArray *group;
	int Dir;
	double Value;
public:
	void Init(int Index,int nGroup, IntArray * group, double Value, int Dir);
};

class Concentrate :
	public Load
{
public:
	Concentrate();
	virtual ~Concentrate();
private:
	IntArray *Nodes;
	int nNode;
	FloatArray *Values;
	int Dir;
public:
	void Init(int Index,int nNode, IntArray * Nodes, FloatArray * Values, int Dir);
};
