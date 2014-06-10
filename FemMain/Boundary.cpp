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
}


Displace::Displace()
{
	nNode = 0;
	Dir = 0;
}


Displace::~Displace()
{
}
