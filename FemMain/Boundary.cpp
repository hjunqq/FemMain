#include "Boundary.h"


Boundary::Boundary()
{
	Index = 0;
}


Boundary::~Boundary()
{
}


Fixed::Fixed()
{
	nNode = 0;
	Dir = 0;
}


Fixed::~Fixed()
{
}


void Fixed::Init(int Index, int nNode, IntArray * Node, int Dir)
{
	this->Index = Index;
	this->nNode = nNode;
	this->Dir = Dir;
	this->Node = new IntArray(*Node);
}


Displace::Displace()
{
	nNode = 0;
	Dir = 0;
}

void Displace::Init(int nNode, int Dir, IntArray *Node, FloatArray *Values)
{
	this->nNode = nNode;
	this->Dir = Dir;
	this->Node = new IntArray(*Node);
	this->Values = new FloatArray(*Values);
}


Displace::~Displace()
{
}

Interact::Interact()
{
	nNode = 0;
	AdjDomain = 0;
}

void Interact::Init(int nNode, int AdjDomain, IntArray *Local, IntArray *Remote)
{
	this->nNode = nNode;
	this->AdjDomain = AdjDomain;
	this->Local = new IntArray(*Local);
	this->Remote = new IntArray(*Remote);
}

Interact::~Interact()
{
}

