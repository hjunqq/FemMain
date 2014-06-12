#pragma once
#include "Array.h"
class Boundary
{
public:
	Boundary();
	virtual ~Boundary();
protected:
	int Index;
};

class Fixed :
	public Boundary
{
public:
	Fixed();
	virtual ~Fixed();
	void Init(int Index, int nNode, IntArray * Node, int Dir);
private:
	int nNode;
	IntArray *Node;
	int Dir;
public:
	// 获得第i个元素
	int at(int i);
	int GetDir();
	int GetnNode();
};

class Displace :
	public Boundary
{
public:
	Displace();
	virtual ~Displace();
	void Init(int nNode, int Dir, IntArray *Node, FloatArray *Values);
private:
	int nNode;
	IntArray *Node;
	FloatArray *Values;
	int Dir;
public:
	int GetnNode();
	int at(int i);
	double atValue(int i);
	int GetDir();
};

class Interact :
	public Boundary
{
public:
	Interact();
	virtual ~Interact();
	void Init(int nNode, int AdjDomain, IntArray *Local, IntArray *Remote);
private:
	int nNode;
	IntArray *Local;
	IntArray *Remote;
	int AdjDomain;
public:
	int GetnNode();
	int GetAdj();
	IntArray * GetLocal();
	IntArray * GetRemote();
};
