#pragma once

#define Elatisic 0

class Material
{
public:
	Material();
	virtual ~Material();
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

