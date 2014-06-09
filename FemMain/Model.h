#pragma once
#include "Array.h"
class Element
{
public:
	Element();
	virtual ~Element();
protected:
	int nNodes;
	int Material;
	int Index;
	int Group;
	IntArray  Node;
	FloatMatrix *Stiff;
	IntArray DegreeOfFreedom;
public:
	// µ¥ÔªÐÎÌ¬
	Element * Type();
	void Init(int nNodes, int Material, int Index, int Group, IntArray Node);
};

