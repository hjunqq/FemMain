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
private:
	int nNode;
	IntArray *Node;
public:
	void Init(int Index, int nNode, IntArray * Node, int Dir);
private:
	int Dir;
};

class Displace :
	public Boundary
{
public:
	Displace();
	virtual ~Displace();
private:
	int nNode;
	IntArray *Node;
	IntArray *Values;
	int Dir;
};

