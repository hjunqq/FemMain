#include "Model.h"


Element::Element()
{
	nNodes = 0;
	Material = 0;
	Index = 0;
	nGaussPoint = 0;
}


Element::~Element()
{
}


 // ��Ԫ��ʼ��
void Element::Init(int nNodes, int Material, int Index, int Group, IntArray *Node)
{
	this->nNodes = nNodes;
	this->Material = Material;
	this->Index = Index;
	this->group = Group;
	this->Nodes = new IntArray(this->nNodes);
	*Nodes = Nodes->Copy(*Node);
	this->DegreeOfFreedom = NULL;
	this->Stiff = NULL;
}


// ����նȾ���
FloatMatrix * Element::BuildStiff()
{
	return NULL;
}


// ��װ��������
FloatMatrix * Element::BuildConstitutiveMatrix()
{
	return NULL;
}


// �����˹��Ӧ��
FloatArray * Element::ComputeStrain(GaussPoint * B)
{
	return NULL;
}


// ����B����
FloatMatrix Element::ComputeBMarix(GaussPoint * B)
{
	return FloatMatrix();
}


// ��ӡ��Ԫ���
void Element::PrintRes()
{
}


// ��õ�Ԫ�ڵ���
int Element::GetnNode()
{
	return 0;
}


// ��õ�Ԫ�ڵ�
IntArray * Element::GetNodeArray()
{
	return Nodes;
}

// ��õ�Ԫ�ڵ�
int Element::GetNode(int i)
{
	return Nodes->at(i);
}


// ��õ�Ԫ��˹����
int Element::GetnGaussPoint()
{
	return nGaussPoint;
}


// ��õ�Ԫ��
int Element::GetGroup()
{
	return group;
}


// ��õ�Ԫ����
int Element::GetMaterial()
{
	return Material;
}


Node::Node()
{
	Index = 0;
}


Node::~Node()
{
}


// �ڵ��ʼ��
void Node::Init(int Index, FloatArray * Coordinate)
{
	this->Index = Index;
	int size;
	size = Coordinate->GetSize();
	Coordinates = new FloatArray(size);
	*Coordinates=Coordinates->Copy(*Coordinate);
}


// ��ӡ�ڵ���Ϣ
void Node::Print()
{
	cout << "Node " << Index ;
	Coordinates->Print();
}


// ��ýڵ�����
double Node::GetCoordinate(int i)
{
	return Coordinates->at(i);
}


// ��ýڵ�����
FloatArray Node::GetCoordinate()
{
	return *Coordinates;
}


// ���λ��
FloatArray & Node::GetDisplacement()
{
	return *Displacement;
	//TODO: insert return statement here
}


// �����Ӧ��
FloatArray & Node::GetPriStess()
{
	return *PrincipleStrain;
	//TODO: insert return statement here
}


// ���Ӧ��
FloatMatrix & Node::GetStress()
{
	return Stress;
	//TODO: insert return statement here
}


// ���Ӧ��
FloatMatrix & Node::GetStrain()
{
	return *Strain;
	//TODO: insert return statement here
}


GaussPoint::GaussPoint()
{
}


GaussPoint::~GaussPoint()
{
}


Group::Group()
{
	Index = 0;
	Material = 0;
	Appear = 0;
	type = 0;
	nElements = 0;
}


Group::~Group()
{
}


// ��ʼ��
void Group::Init(int Index, int nElements, int type , int Mat, int Dof)
{
	this->Index = Index;
	this->nElements = nElements;
	this->type = type;
	this->Material = Mat;
	this->nDof = Dof;
	Elements = new IntArray(nElements);
}


// ��䵥Ԫ
void Group::FillElement(IntArray * ElementList)
{
	*Elements =Elements->Copy(* ElementList);
}


// ��ȡ��i����Ԫ��
int Group::GetElement(int i)
{
	return Elements->at(i);
}




// ��ȡ��Ԫ����
int Group::GetnElements()
{
	return nElements;
}

int Group::GetType()
{
	return type;
}
// ��ȡ���Ϻ�
int Group::GetMaterial()
{
	return Material;
}


// �����Ƿ����
bool & Group::IsAppear()
{
	return Appear;
	//TODO: insert return statement here
}


Quadr::Quadr()
{
	this->type = Quadrilateral;
	this->nNodes = 4;
	this->Nodes = new IntArray(4);
}


Quadr::~Quadr()
{
}


void Element::Print()
{
}


void Quadr::Print()
{
	cout << "Quadr Element "<<setw(5)<<this->Index;
	for (int i = 0; i < 4; i++)
	{
		cout <<setw(5)<< this->Nodes->at(i);
	}
	cout << endl;
}

Line::Line()
{
	this->type = Linear;
	this->nNodes = 2;
	this->Nodes = new IntArray(2);
}
Line::~Line()
{
}
Line::Line(const Line & L)
{
	nNodes = L.nNodes;
	Material = L.Material;
	Index = L.Index;
	group = L.group;
	type = L.type;
	Nodes = new IntArray(nNodes);
	Nodes = L.Nodes;
}


void Line::Print()
{
	cout << "I am a Line Element" << endl;
}

int & Line::AtAdjElem()
{
	return AdjElem;
}