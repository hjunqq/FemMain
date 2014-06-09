#include "Model.h"


Element::Element()
{
	nNodes = 0;
	Material = 0;
	Index = 0;
}


Element::~Element()
{
}


// 单元形态
Element * Element::Type()
{
	return NULL;
}
 // 单元初始化
void Element::Init(int nNodes, int Material, int Index, int Group, IntArray Node)
{
	this->nNodes = nNodes;
	this->Material = Material;
	this->Index = Index;
	this->Group = Group;
	this->Node = Node;
	this->DegreeOfFreedom = NULL;
	this->Stiff = NULL;
}
