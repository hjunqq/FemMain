#pragma once
#include "Array.h"
#include "Material.h"

#define Linear			2
#define Quadrilateral	4
const double Gauss3[] = { -0.774596669241483, 0.0, 0.774596669241483 };
const double Weight3[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
const double Gauss2[] = { -0.5773502692, 0.5773502692 };
const double Weight2[] = { 1.0, 1.0};
const double P4ksicor[] = { -0.577350269, 0.577350269, 0.577350269, -0.577350269 };
const double P4etacor[] = { -0.577350269, -0.577350269, 0.577350269, 0.577350269 };

class Node
{
public:
	Node();
	virtual ~Node();
private:
	int Index;
	int Count;
	FloatArray Coordinates;
	FloatArray Displacement;
	FloatArray PrincipleStrain;
	FloatArray  Stress;
	FloatArray  Strain;
public:
	// 节点初始化
	void Init(int Index, FloatArray & Coordinate);
	// 打印节点信息
	void Print();
	void ResetCount();
	void AddCount();
	// 获得节点坐标
	double GetCoordinate(int i);
	// 获得节点坐标
	FloatArray  GetCoordinate();
	// 获得位移
	FloatArray  GetDisplacement();
	// 获得主应力
	FloatArray  GetPriStess();
	// 获得应力
	FloatArray  GetStress();
	// 获得应变
	FloatArray  GetStrain();

	void SetDisplacement(FloatArray &NodeDisplacement);
	void SetStress(FloatArray & NodeStress);
	void SetStrain(FloatArray & NodeStrain);
	Node & operator =(const Node & N);
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
	int mat;
	int Index;
	int group;
	int Dof;
	int type;
	double Det;
	Material Mat;
	IntArray  Nodes;
	FloatMatrix Coors;
	FloatArray Shape;
	FloatMatrix DShape;
	FloatMatrix DShapeX;
	FloatMatrix Jacobi;
	FloatMatrix InvJacobi;
	FloatMatrix Stiff,Mass,Damp;
	FloatArray Strain;
	FloatArray Stress;
	FloatMatrix BMatrix;
	FloatMatrix BMatrixBig;
	FloatMatrix DMatrix;
	IntArray DegreeOfFreedom;
	FloatArray Displacement;
	FloatMatrix ConstitutiveMatrix;
	int nGaussPoint;
	GaussPoint *GaussPointArray = NULL;
	FloatArray *GaussStrain = NULL;
	FloatArray *GaussStress = NULL;
	FloatArray *NodeStrain = NULL;
	FloatArray *NodeStress = NULL;
public:
	void Init(int nNodes, int mat, int Index, int Group, int Dof, IntArray *Node);
	// 计算J
	virtual FloatMatrix  ComputeJacobi(GaussPoint * B);
	// 计算刚度矩阵
	virtual FloatMatrix  ComputeStiff();
	// 组装本构矩阵
	virtual FloatMatrix  ComputeConstitutiveMatrix();
	// 计算高斯点应变
	virtual FloatArray  ComputeStress();
	// 计算B矩阵
	virtual FloatMatrix  ComputeBMarix(GaussPoint * B);
	// 打印单元结果
	void PrintRes();
	// 获得单元节点数
	int GetnNode();
	// 获得单元节点
	IntArray  GetNodeArray();
	// 获得单元节点
	int GetNode(int i);
	// 获得单元高斯点数
	int GetnGaussPoint();
	// 获得单元组
	int GetGroup();
	// 获得单元材料
	int GetMaterial();
	void SetMaterial(Material & Mat);

	void SetCoor(FloatMatrix & Coor);
	void SetResult(FloatArray & Result);
	void SetInitialDisplacement(FloatArray &InitialDisplacement);

	int GetIndex();
	void FillDof(IntArray & DegreeOfFreedom);
	IntArray GetDof();
	FloatMatrix GetStiff();
	FloatMatrix GetMass();
	FloatArray GetShape();
	FloatArray GetStrain(int inode);
	FloatArray GetStress(int inode);
	FloatArray GetDisplacement(int inode);
	double GetDet();
	virtual void Print();
};

class Quadr :
	public Element
{
public:
	Quadr();
	virtual ~Quadr();
	void Print();
	Quadr & operator=(const Quadr &Q);
	FloatMatrix  ComputeStiff();
	FloatMatrix  ComputeJacobi(GaussPoint & B);
	FloatMatrix  ComputeBMarix(int inode);
	FloatMatrix  ComputeBMatrix();
	FloatMatrix  ComputeConstitutiveMatrix();
	FloatArray  ComputeStress();
	FloatMatrix ComputeMassMatrix();
};
class Line :
	public Element
{
private:
	int AdjElem;
public:
	Line();
	Line(const Line & L);
	virtual ~Line();
	void Print();
	int & AtAdjElem();
	FloatArray ComputeShape(GaussPoint *B);
	Line & operator=(const Line &L);
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
	int nDof;
	double Alpha, Beta;
	Quadr ** Quadrs; 
public:
	// 初始化
	void Init(int Index, int nElements, int type,int Mat,int Dof);
	// 填充单元
	void FillElement(IntArray * ElementList);
	void FillElement(Quadr * Quadrs);
	void FillDampPara(double Alpha, double Beta);
	// 获取第i个单元号
	Quadr ** GetElement(Quadr &);
	// 获取单元个数
	int GetnElements();
	// 获取材料号
	int GetMaterial();
	int GetType();
	int GetDof();
	double GetDampAlpha();
	double GetDampBeta();
	// 设置是否出现
	bool  IsAppear();
};



