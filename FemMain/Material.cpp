#include "Material.h"


Material::Material()
{
	Index = 0;
	type = 0;
	Young = 0.0;
	Possion = 0.0;
	Density = 0.0;
}

Material::Material(const Material &Mat)
{
	Index = Mat.Index;
	type = Mat.type;
	Density = Mat.Density;
	Possion = Mat.Possion;
	Young = Mat.Young;
}

Material::~Material()
{
}


double Material::GetYoung()
{
	return Young;
}


double Material::GetPossion()
{
	return Possion;
}


double Material::GetDensity()
{
	return Density;
}


// ³õÊ¼»¯
void Material::Init(int Index, int type, double Young, double Possion, double Density)
{
	this->Index = Index;
	this->type = type;
	this->Young = Young;
	this->Possion = Possion;
	this->Density = Density;
}
