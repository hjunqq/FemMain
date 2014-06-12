#include "Model.h"
#include "FemMain.h"
extern FemMain Fem1;

Element::Element()
{
	nNodes = 0;
	mat = 0;
	Index = 0;
	nGaussPoint = 0;
}


Element::~Element()
{
}


 // ��Ԫ��ʼ��
void Element::Init(int nNodes, int Material, int Index, int Group,int Dof, IntArray *Node)
{
	this->nNodes = nNodes;
	this->mat = Material;
	this->Index = Index;
	this->group = Group;
	this->Dof = Dof;
	this->Nodes = new IntArray(*Node);
	this->DegreeOfFreedom = NULL;
	this->Stiff = NULL;
}

FloatMatrix *Element::ComputeJacobi(GaussPoint *B)
{
	return NULL;
}

FloatMatrix *Element::ComputeStiff()
{
	return NULL;
}

FloatMatrix *Element::ComputeConstitutiveMatrix()
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
	return mat;
}

void Element::SetMaterial(Material *Mat)
{
	this->Mat =new Material( *Mat);
}

void Element::SetCoor(FloatArray **Coor)
{
	this->Coors = new FloatArray*[nNodes];
	for (int i = 0; i < nNodes; i++)
	{
		this->Coors[i] = new FloatArray(*Coor[i]);
	}
}

void Element::Print()
{
}

int Element::GetIndex()
{
	return Index;
}

void Element::FillDof(IntArray * DegreeOfFreedom)
{
	int Node;
	this->DegreeOfFreedom = new IntArray(nNodes*Dof);
	for (int inode = 0; inode < nNodes; inode++)
	{
		Node = Nodes->at(inode);
		for (int idof = 0; idof < Dof; idof++)
		{
			this->DegreeOfFreedom->at(inode*Dof + idof) = DegreeOfFreedom->at(Node*Dof + idof);
		}
	}
}

IntArray *Element::GetDof()
{
	return DegreeOfFreedom;
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
	Coordinates = new FloatArray(*Coordinate);
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

Node & Node::operator=(const Node & N)
{
	Node T;
	T.Index = N.Index;
	T.Coordinates = N.Coordinates;
	T.Displacement = N.Displacement;
	T.Stress = N.Stress;
	T.Strain = N.Strain;
	T.PrincipleStrain = N.PrincipleStrain;
	return T;
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
	Quadrs = new Quadr *[nElements];
}


// ��䵥Ԫ
void Group::FillElement(IntArray * ElementList)
{
	*Elements =Elements->Copy(* ElementList);
}

void Group::FillElement(Quadr *Q)
{
	for (int ielem = 0; ielem < nElements;ielem++)
	{
		Quadrs[ielem] = &Q[ielem];
	}
}

// ��ȡ��i����Ԫ
Quadr * Group::GetElement(Quadr & )
{
	return *Quadrs;
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

int Group::GetDof()
{
	return nDof;
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

FloatMatrix * Quadr::ComputeJacobi(GaussPoint * B)
{
	FloatArray **Coor;
	double Young, Possion;
	FloatMatrix *DN,*DX;
	return NULL;
}
FloatMatrix * Quadr::BuildStiff()
{
	return NULL;
}

Quadr::~Quadr()
{
}


Quadr & Quadr::operator=(const Quadr &Q)
{
	Quadr T;
	T.Index = Q.Index;
	T.ConstitutiveMatrix = Q.ConstitutiveMatrix;
	T.DegreeOfFreedom = Q.DegreeOfFreedom;
	T.GaussPointArray = Q.GaussPointArray;
	T.DegreeOfFreedom = Q.DegreeOfFreedom;
	T.group = Q.group;
	T.mat = Q.mat;
	T.nGaussPoint = Q.nGaussPoint;
	T.nNodes = Q.nNodes;
	T.Nodes = Q.Nodes;
	T.Stiff = Q.Stiff;
	T.type = Q.type;
	return T;
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
	mat = L.mat;
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

Line & Line::operator=(const Line &L)
{
	Line T;
	T.Index = L.Index;
	T.ConstitutiveMatrix = L.ConstitutiveMatrix;
	T.DegreeOfFreedom = L.DegreeOfFreedom;
	T.GaussPointArray = L.GaussPointArray;
	T.DegreeOfFreedom = L.DegreeOfFreedom;
	T.group = L.group;
	T.mat = L.mat;
	T.nGaussPoint = L.nGaussPoint;
	T.nNodes = L.nNodes;
	T.Nodes = L.Nodes;
	T.Stiff = L.Stiff;
	T.type = L.type;
	T.AdjElem = L.AdjElem;
	return T;
}