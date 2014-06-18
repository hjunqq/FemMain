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
	int nNode;
	double StartV;
	double EndV;
	int Dir;
	double StartC;
	double EndC;
	Line **Lines;
	Quadr **Quadrs;
public:
	void Init(int Index,int nEle,int nNode, double StartC, double StartV, double EndC, double EndV, int Dir);
	void Set(Line *Lines);
	void Set(Quadr *Quadrs);
	int GetnElem();
	int GetnNode();
	int GetDir();
	double GetStartV();
	double GetEndV();
	double GetStartC();
	double GetEndC();
	Line *GetLines();
	Quadr *GetQuadr();
	
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
	IntArray group;
	int Dir;
	double Value;
public:
	void Init(int Index,int nGroup, IntArray & group, double Value, int Dir);
	int GetnGroup();
	int GetDir();
	double GetValue();
	IntArray GetGroup();
};

class Concentrate :
	public Load
{
public:
	Concentrate();
	virtual ~Concentrate();
private:
	IntArray Nodes;
	int nNode;
	FloatArray Values;
	int Dir;
public:
	void Init(int Index,int nNode, IntArray & Nodes, FloatArray & Values, int Dir);
	int GetnNode();
	int GetDir();
	IntArray GetNodes();
	FloatArray GetValues();
};
