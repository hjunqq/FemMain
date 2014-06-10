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
};

class Face :
	public Load
{
public:
	Face();
	virtual ~Face();
private:
	int nElem;
	Element *Elems;
	double StratV;
	double EndV;
	int Dir;
	double StartC;
	double EndC;
public:
	void Init(int Index, double StartC, double StartV, double EndC, double EndV, int Dir);
};

class Volumn :
	public Load
{
public:
	Volumn();
	virtual ~Volumn();
private:
	int Index;
	int group;
	int Dir;
	double Value;
public:
	void Init(int Index, int group, double Value, int Dir);
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
	void Init(int Index, IntArray * Nodes, FloatArray * Values, int Dir);
};
