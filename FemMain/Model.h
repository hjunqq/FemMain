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
	// �ڵ��ʼ��
	void Init(int Index, FloatArray * Coordinate);
	// ��ӡ�ڵ���Ϣ
	void Print();
	// ��ýڵ�����
	double GetCoordinate(int i);
	// ��ýڵ�����
	FloatArray GetCoordinate();
	// ���λ��
	FloatArray & GetDisplacement();
	// �����Ӧ��
	FloatArray & GetPriStess();
	// ���Ӧ��
	FloatMatrix & GetStress();
	// ���Ӧ��
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
	// ����նȾ���
	virtual FloatMatrix * BuildStiff();
	// ��װ��������
	virtual FloatMatrix * BuildConstitutiveMatrix();
	// �����˹��Ӧ��
	virtual FloatArray * ComputeStrain(GaussPoint * B);
	// ����B����
	virtual FloatMatrix ComputeBMarix(GaussPoint * B);
	// ��ӡ��Ԫ���
	void PrintRes();
	// ��õ�Ԫ�ڵ���
	int GetnNode();
	// ��õ�Ԫ�ڵ�
	IntArray GetNodeArray();
	// ��õ�Ԫ�ڵ�
	int GetNode(int i);
	// ��õ�Ԫ��˹����
	int GetnGaussPoint();
	// ��õ�Ԫ��
	int GetGroup();
	// ��õ�Ԫ����
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
	// ��ʼ��
	void Init(int Index, int nElements, int type);
	// ��䵥Ԫ
	void FillElement(IntArray * ElementList);
	// ��ȡ��i����Ԫ��
	int GetElement(int i);
	// ��ȡ��Ԫ����
	int GetnElements();
	// ��ȡ���Ϻ�
	int GetMaterial();
	// �����Ƿ����
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

