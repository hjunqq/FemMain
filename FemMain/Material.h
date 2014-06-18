#pragma once

#define Elastic 0

class Material
{
public:
	Material();
	virtual ~Material();
	Material(const Material & Mat);
protected:
	int Index;
	int type;
	double Young;
	double Possion;
	double Density;
public:
	double GetYoung();
	double GetPossion();
	double GetDensity();
	// ≥ı ºªØ
	virtual void Init(int Index, int type, double Young, double Possion, double Density);
};

