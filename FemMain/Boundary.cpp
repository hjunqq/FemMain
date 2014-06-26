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


void Fixed::Init(int Index, int nNode, IntArray & Node, int Dir)
{
	this->Index = Index;
	this->nNode = nNode;
	this->Dir = Dir;
	this->Node = Node;
}


Displace::Displace()
{
	nNode = 0;
	Dir = 0;
}

void Displace::Init(int nNode, int Dir, IntArray &Node, FloatArray &Values)
{
	this->nNode = nNode;
	this->Dir = Dir;
	this->Node =Node;
	this->Values =Values;
}


Displace::~Displace()
{
}

Interact::Interact()
{
	nNode = 0;
	AdjDomain = 0;
}

void Interact::Init(int nNode, int AdjDomain, IntArray &Local, IntArray &Remote)
{
	this->nNode = nNode;
	this->AdjDomain = AdjDomain;
	this->Local = Local;
	this->Remote = Remote;
}

Interact::~Interact()
{
}



// 获得第i个元素
int Fixed::at(int i)
{
	return Node.at(i);
}


int Fixed::GetDir()
{
	return Dir;
}


int Fixed::GetnNode()
{
	return nNode;
}


int Displace::GetnNode()
{
	return nNode;
}


int Displace::at(int i)
{
	return Node.at(i);
}


double Displace::atValue(int i)
{
	return Values.at(i);
}

IntArray Displace::GetNode()
{
	return Node;
}

FloatArray Displace::GetValue()
{
	return Values;
}


int Displace::GetDir()
{
	return Dir;
}


int Interact::GetnNode()
{
	return nNode;
}


int Interact::GetAdj()
{
	return AdjDomain;
}


IntArray  Interact::GetLocal()
{
	return Local;
}


IntArray  Interact::GetRemote()
{
	return Remote;
}
