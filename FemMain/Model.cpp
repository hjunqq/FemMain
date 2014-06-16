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


 // 单元初始化
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

// 计算高斯点应变
FloatArray * Element::ComputeStrain(GaussPoint * B)
{
	return NULL;
}


// 计算B矩阵
FloatMatrix ** Element::ComputeBMarix(GaussPoint * B)
{
	return NULL;
}



// 打印单元结果
void Element::PrintRes()
{
}


// 获得单元节点数
int Element::GetnNode()
{
	return 0;
}


// 获得单元节点
IntArray * Element::GetNodeArray()
{
	return Nodes;
}

// 获得单元节点
int Element::GetNode(int i)
{
	return Nodes->at(i);
}


// 获得单元高斯点数
int Element::GetnGaussPoint()
{
	return nGaussPoint;
}


// 获得单元组
int Element::GetGroup()
{
	return group;
}


// 获得单元材料
int Element::GetMaterial()
{
	return mat;
}

void Element::SetMaterial(Material *Mat)
{
	this->Mat =new Material( *Mat);
}

void Element::SetCoor(FloatMatrix *Coor)
{
	this->Coors = new FloatMatrix(*Coor);
}

void Element::SetResult(FloatArray *Result)
{
	if (Displacement == NULL)
	{
		Displacement = new FloatArray(Dof);
	}
	else
	{
		Displacement->Clear();
	}
	for (int iDof = 0; iDof < Dof; iDof++)
	{
		int DofIdx = DegreeOfFreedom->at(iDof);
		if (DofIdx != 0)
		{
			Displacement->at(iDof - 1) = Result->at(DofIdx - 1);
		}
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

FloatMatrix *Element::GetStiff()
{
	return Stiff;
}
FloatArray *Element::GetShape()
{
	return Shape;
}
double Element::GetDet()
{
	return Det;
}
Node::Node()
{
	Index = 0;
}


Node::~Node()
{
}


// 节点初始化
void Node::Init(int Index, FloatArray * Coordinate)
{
	this->Index = Index;
	int size;
	size = Coordinate->GetSize();
	Coordinates = new FloatArray(*Coordinate);
}


// 打印节点信息
void Node::Print()
{
	cout << "Node " << Index ;
	Coordinates->Print();
}


// 获得节点坐标
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

// 获得节点坐标
FloatArray & Node::GetCoordinate()
{
	return *Coordinates;
}


// 获得位移
FloatArray & Node::GetDisplacement()
{
	return *Displacement;
	//TODO: insert return statement here
}


// 获得主应力
FloatArray & Node::GetPriStess()
{
	return *PrincipleStrain;
	//TODO: insert return statement here
}


// 获得应力
FloatMatrix & Node::GetStress()
{
	return *Stress;
	//TODO: insert return statement here
}


// 获得应变
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


// 初始化
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


// 填充单元
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

// 获取第i个单元
Quadr * Group::GetElement(Quadr & )
{
	return *Quadrs;
}




// 获取单元个数
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

// 获取材料号
int Group::GetMaterial()
{
	return Material;
}


// 设置是否出现
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
	if (this->Coors == NULL)
	{
		cout << "No Coordiantes" << endl;
	}
	Shape = new FloatArray(4);
	double ksi, eta;
	ksi = B->GetCoordinate(0);
	eta = B->GetCoordinate(1);
	Shape->at(0) = 0.25*(1 + ksi)*(1 + eta);
	Shape->at(1) = 0.25*(1 - ksi)*(1 + eta);
	Shape->at(2) = 0.25*(1 - ksi)*(1 - eta);
	Shape->at(3) = 0.25*(1 + ksi)*(1 - eta);
	DShape = new FloatMatrix(2, 4);
	DShape->at(0, 0) = 0.25*(1 + eta);
	DShape->at(1, 0) = 0.25*(1 + ksi);
	DShape->at(0, 1) = -0.25*(1 + eta);
	DShape->at(1, 1) = 0.25*(1 - ksi);
	DShape->at(0, 2) = -0.25*(1 - eta);
	DShape->at(1, 2) = -0.25*(1 - ksi);
	DShape->at(0, 3) = 0.25*(1 - eta);
	DShape->at(1, 3) = -0.25*(1 + ksi);
	Jacobi = new FloatMatrix(2, 2);
	//DShape->Print();
	//Coors->Print();
	Jacobi =& DShape->Mult(*Coors);
	//Jacobi->Print();
	Det = Jacobi->Determinant();
	InvJacobi = new FloatMatrix(2, 2);
	InvJacobi = &Jacobi->Inverse();
	//InvJacobi->Print();
	DShapeX = new FloatMatrix(2, 4);
	DShapeX = &InvJacobi->Mult(*DShape);
	return NULL;
}
FloatMatrix ** Quadr::ComputeBMarix()
{
	BMatrix = new FloatMatrix*[4];
	for (int inode = 0; inode < 4; inode++)
	{
		BMatrix[inode] = new FloatMatrix(3, 2);
		BMatrix[inode]->at(0, 0) = DShape->at(0, inode);
		BMatrix[inode]->at(1, 1) = DShape->at(1, inode);
		BMatrix[inode]->at(2, 0) = DShape->at(1, inode);
		BMatrix[inode]->at(2, 1) = DShape->at(0, inode);
	}
	return BMatrix;
}
FloatMatrix * Quadr::ComputeConstitutiveMatrix()
{
	FloatMatrix *T;
	T = new FloatMatrix(3, 3);
	double Young, Possion;
	Young = Mat->GetYoung();
	Possion = Mat->GetPossion();
	double D1, D2, D3;
	D1 = (1 - Possion)*Young / ((1 + Possion)*(1 - 2 * Possion));
	D2 = Possion*Young / ((1 + Possion)*(1 - 2 * Possion));
	D3 = Young / (2 * (1 - Possion));
	T->at(0, 0) = D1;
	T->at(1, 0) = D2;
	T->at(0, 1) = D2;
	T->at(1, 1) = D1;
	T->at(2, 2) = D3;
	return T;
}
FloatMatrix * Quadr::ComputeStiff()
{
	FloatArray *GaussCoor;
	GaussCoor = new FloatArray(2);
	GaussPoint *B = new GaussPoint;
	DMatrix = new FloatMatrix(3, 3);
	FloatMatrix *G = new FloatMatrix(2, 2);
	FloatMatrix *BT=new FloatMatrix(3, 2);
	FloatMatrix *DT = new FloatMatrix(3, 3);
	Stiff = new FloatMatrix(8, 8);
	double W1,W2;
	for (int i = 0; i < 3; i++)
	{
		GaussCoor->at(0) = Gauss3[i];
		W1 = Weight3[i];
		for (int j = 0; j < 3; j++)
		{
			GaussCoor->at(1) = Gauss3[j];
			B->Init(0, GaussCoor);
			ComputeJacobi(B);
			BMatrix = ComputeBMarix();
			DMatrix = ComputeConstitutiveMatrix();
			for (int inode = 0; inode < 4; inode++)
			{
				for (int jnode = 0; jnode < 4; jnode++)
				{
					BT = &BMatrix[inode]->Trans();
					DT = &BT->Mult(*DMatrix);
					G = &DT->Mult(*BMatrix[jnode]);
					W2 = Weight3[j];
					W2 = W1*W2;
					G = &G->Mult(W2);
					Stiff->at(2*inode, 2*jnode) += G->at(0, 0);
					Stiff->at(2*inode+1, 2*jnode) += G->at(1, 0);
					Stiff->at(2*inode, 2*jnode+1) += G->at(0, 1);
					Stiff->at(2*inode+1, 2*jnode+1) += G->at(1, 1);
				}
			}
		}
	}
	//Stiff->Print();
	return Stiff;
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

FloatArray *Line::ComputeShape(GaussPoint *B)
{
	Shape = new FloatArray(2);
	double ksi;
	ksi = B->GetCoordinate(0);
	Shape->at(0) = 0.5*(1 - ksi);
	Shape->at(1) = 0.5*(1 + ksi);
	return Shape;
}