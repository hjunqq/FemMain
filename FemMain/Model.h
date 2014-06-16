#pragma once
#include "Array.h"
#include "Material.h"

#define Linear			2
#define Quadrilateral	4
const double Gauss3[] = { -0.774596669241483, 0.0, 0.774596669241483 };
const double Weight3[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
const double Gauss2[] = { -0.5773502692, 0.5773502692 };
const double Weight2[] = { 1.0, 1.0};

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
	FloatMatrix * Stress;
	FloatMatrix *Strain;
public:
	// �ڵ��ʼ��
	void Init(int Index, FloatArray * Coordinate);
	// ��ӡ�ڵ���Ϣ
	void Print();
	// ��ýڵ�����
	double GetCoordinate(int i);
	// ��ýڵ�����
	FloatArray & GetCoordinate();
	// ���λ��
	FloatArray & GetDisplacement();
	// �����Ӧ��
	FloatArray & GetPriStess();
	// ���Ӧ��
	FloatMatrix & GetStress();
	// ���Ӧ��
	FloatMatrix & GetStrain();
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
	Material *Mat;
	IntArray  *Nodes;
	FloatMatrix *Coors;
	FloatArray *Shape;
	FloatMatrix *DShape;
	FloatMatrix *DShapeX;
	FloatMatrix *Jacobi;
	FloatMatrix *InvJacobi;
	FloatMatrix *Stiff;
	FloatMatrix **BMatrix;
	FloatMatrix *DMatrix;
	IntArray *DegreeOfFreedom;
	FloatArray *Displacement;
	FloatMatrix *ConstitutiveMatrix;
	int nGaussPoint;
	GaussPoint **GaussPointArray;
public:
	void Init(int nNodes, int mat, int Index, int Group, int Dof, IntArray *Node);
	// ����J
	virtual FloatMatrix * ComputeJacobi(GaussPoint * B);
	// ����նȾ���
	virtual FloatMatrix * ComputeStiff();
	// ��װ��������
	virtual FloatMatrix * ComputeConstitutiveMatrix();
	// �����˹��Ӧ��
	virtual FloatArray * ComputeStrain(GaussPoint * B);
	// ����B����
	virtual FloatMatrix ** ComputeBMarix(GaussPoint * B);
	// ��ӡ��Ԫ���
	void PrintRes();
	// ��õ�Ԫ�ڵ���
	int GetnNode();
	// ��õ�Ԫ�ڵ�
	IntArray * GetNodeArray();
	// ��õ�Ԫ�ڵ�
	int GetNode(int i);
	// ��õ�Ԫ��˹����
	int GetnGaussPoint();
	// ��õ�Ԫ��
	int GetGroup();
	// ��õ�Ԫ����
	int GetMaterial();
	void SetMaterial(Material *Mat);

	void SetCoor(FloatMatrix *Coor);
	void SetResult(FloatArray *Result);

	int GetIndex();
	void FillDof(IntArray * DegreeOfFreedom);
	IntArray * GetDof();
	FloatMatrix *GetStiff();
	FloatArray *GetShape();
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
	FloatMatrix * ComputeStiff();
	FloatMatrix * ComputeJacobi(GaussPoint * B);
	FloatMatrix ** ComputeBMarix();
	FloatMatrix * ComputeConstitutiveMatrix();
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
	FloatArray *ComputeShape(GaussPoint *B);
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
	Quadr ** Quadrs; 
public:
	// ��ʼ��
	void Init(int Index, int nElements, int type,int Mat,int Dof);
	// ��䵥Ԫ
	void FillElement(IntArray * ElementList);
	void FillElement(Quadr * Quadrs);
	// ��ȡ��i����Ԫ��
	Quadr * GetElement(Quadr &);
	// ��ȡ��Ԫ����
	int GetnElements();
	// ��ȡ���Ϻ�
	int GetMaterial();
	int GetType();
	int GetDof();
	// �����Ƿ����
	bool & IsAppear();
};



