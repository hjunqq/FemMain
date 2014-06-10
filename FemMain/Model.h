#pragma once
#include "Array.h"

#define Linear			2
#define Quadrilateral	4

class Node
{
public:
	Node();
	virtual ~Node();
private:
	int Index;
	FloatArray *Coordinates;
	FloatArray *Displacement;
	FloatArray *PrincipleStrain;
	FloatMatrix Stress;
	FloatMatrix *Strain;
public:
	// 节点初始化
	void Init(int Index, FloatArray * Coordinate);
	// 打印节点信息
	void Print();
	// 获得节点坐标
	double GetCoordinate(int i);
	// 获得节点坐标
	FloatArray GetCoordinate();
	// 获得位移
	FloatArray & GetDisplacement();
	// 获得主应力
	FloatArray & GetPriStess();
	// 获得应力
	FloatMatrix & GetStress();
	// 获得应变
	FloatMatrix & GetStrain();
};

class GaussPoint :
	public Node
{
public:
	GaussPoint();
	virtual ~GaussPoint();
};

class Element
{
public:
	Element();
	virtual ~Element();
protected:
	int nNodes;
	int Material;
	int Index;
	int group;
	int type;
	IntArray  Nodes;
	FloatMatrix *Stiff;
	IntArray DegreeOfFreedom;
	FloatMatrix *ConstitutiveMatrix;
	int nGaussPoint;
	GaussPoint **GaussPointArray;
public:
	void Init(int nNodes, int Material, int Index, int Group, IntArray Node);
	// 计算刚度矩阵
	virtual FloatMatrix * BuildStiff();
	// 组装本构矩阵
	virtual FloatMatrix * BuildConstitutiveMatrix();
	// 计算高斯点应变
	virtual FloatArray * ComputeStrain(GaussPoint * B);
	// 计算B矩阵
	virtual FloatMatrix ComputeBMarix(GaussPoint * B);
	// 打印单元结果
	void PrintRes();
	// 获得单元节点数
	int GetnNode();
	// 获得单元节点
	IntArray GetNodeArray();
	// 获得单元节点
	int GetNode(int i);
	// 获得单元高斯点数
	int GetnGaussPoint();
	// 获得单元组
	int GetGroup();
	// 获得单元材料
	int GetMaterial();
	virtual void Print();
};

class Group
{
public:
	Group();
	virtual ~Group();
private:
	int Index;
	int Material;
	bool Appear;
	int type;
	IntArray *Elements;
	int nElements;
public:
	// 初始化
	void Init(int Index, int nElements, int type);
	// 填充单元
	void FillElement(IntArray * ElementList);
	// 获取第i个单元号
	int GetElement(int i);
	// 获取单元个数
	int GetnElements();
	// 获取材料号
	int GetMaterial();
	// 设置是否出现
	bool & IsAppear();
};

class Quadr :
	public Element
{
public:
	Quadr();
	virtual ~Quadr();
	void Print();
};
class Line :
	public Element
{
public:
	Line();
	virtual ~Line();
	void Print();
};

