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
	this->Nodes = *Node;
}

FloatMatrix Element::ComputeJacobi(GaussPoint *B)
{
	return Jacobi;
}

FloatMatrix Element::ComputeStiff()
{
	return Stiff;
}

FloatMatrix Element::ComputeConstitutiveMatrix()
{
	return ConstitutiveMatrix;
}

// 计算高斯点应变
FloatArray  Element::ComputeStress()
{
	return Strain;
}


// 计算B矩阵
FloatMatrix  Element::ComputeBMarix(GaussPoint * B)
{
	return BMatrix;
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
IntArray  Element::GetNodeArray()
{
	return Nodes;
}

// 获得单元节点
int Element::GetNode(int i)
{
	return Nodes.at(i);
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

void Element::SetMaterial(Material & Mat)
{
	this->Mat =Mat;
}

void Element::SetCoor(FloatMatrix &Coor)
{
	this->Coors= Coor;
}

void Element::SetResult(FloatArray & Result)
{
	for (int iDof = 0; iDof < Dof*nNodes; iDof++)
	{
		int DofIdx = DegreeOfFreedom.at(iDof);
		if (DofIdx != 0)
		{
			Displacement.at(iDof) = Result.at(DofIdx - 1);
		}
	}
}

FloatArray Element::GetStrain(int inode)
{
	return NodeStrain[inode];
}
FloatArray Element::GetStress(int inode)
{
	return NodeStress[inode];
}

void Element::Print()
{
}

int Element::GetIndex()
{
	return Index;
}

void Element::FillDof(IntArray & DegreeOfFreedom)
{
	int Node;
	this->DegreeOfFreedom.SetSize(nNodes*Dof);
	for (int inode = 0; inode < nNodes; inode++)
	{
		Node = Nodes.at(inode);
		for (int idof = 0; idof < Dof; idof++)
		{
			this->DegreeOfFreedom.at(inode*Dof + idof) = DegreeOfFreedom.at(Node*Dof + idof);
		}
	}
}

IntArray Element::GetDof()
{
	return DegreeOfFreedom;
}

FloatMatrix Element::GetStiff()
{
	return Stiff;
}
FloatArray Element::GetShape()
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
	Count = 0;
}


Node::~Node()
{
}


// 节点初始化
void Node::Init(int Index, FloatArray & Coordinate)
{
	this->Index = Index;
	int size;
	size = Coordinate.GetSize();
	Coordinates = Coordinate;
}


// 打印节点信息
void Node::Print()
{
	cout << "Node " << Index ;
	Coordinates.Print();
}

void Node::ResetCount()
{
	Count = 0;
}
void Node::AddCount()
{
	Count++;
}
// 获得节点坐标
double Node::GetCoordinate(int i)
{
	return Coordinates.at(i);
}

Node & Node::operator=(const Node & N)
{
	Index = N.Index;
	Coordinates = N.Coordinates;
	Displacement = N.Displacement;
	Stress = N.Stress;
	Strain = N.Strain;
	PrincipleStrain = N.PrincipleStrain;
	return *this;
}

// 获得节点坐标
FloatArray Node::GetCoordinate()
{
	return Coordinates;
}


// 获得位移
FloatArray Node::GetDisplacement()
{
	return Displacement;
	//TODO: insert return statement here
}


// 获得主应力
FloatArray Node::GetPriStess()
{
	return PrincipleStrain;
	//TODO: insert return statement here
}


// 获得应力
FloatMatrix Node::GetStress()
{
	return Stress;
	//TODO: insert return statement here
}


// 获得应变
FloatMatrix Node::GetStrain()
{
	return Strain;
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
Quadr ** Group::GetElement(Quadr & )
{
	return Quadrs;
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
bool Group::IsAppear()
{
	return Appear;
	//TODO: insert return statement here
}


Quadr::Quadr()
{
	this->type = Quadrilateral;
	this->nNodes = 4;
	this->Nodes.SetSize(4);
	nGaussPoint = 4;
	Shape.SetSize(4);
	DShape.SetSize(2, 4);
	Jacobi.SetSize(2, 2);
	DShapeX.SetSize(2, 4);
	BMatrix.SetSize(3, 2);
	BMatrixBig.SetSize(3, 8);
	DMatrix.SetSize(3, 3);
	Stiff.SetSize(8, 8);
	Displacement.SetSize(8);
	GaussStrain = new FloatArray[nGaussPoint];
	GaussStress = new FloatArray[nGaussPoint];
	NodeStrain = new FloatArray[nNodes];
	NodeStress = new FloatArray[nNodes];
	for (int i = 0; i < 4; i++)
	{
		GaussStrain[i].SetSize(3);
		GaussStress[i].SetSize(3);
		NodeStrain[i].SetSize(3);
		NodeStress[i].SetSize(3);
	}
}

FloatMatrix Quadr::ComputeJacobi(GaussPoint & B)
{
	double ksi, eta;
	ksi = B.GetCoordinate(0);
	eta = B.GetCoordinate(1);
	Shape.at(0) = 0.25*(1 + ksi)*(1 + eta);
	Shape.at(1) = 0.25*(1 - ksi)*(1 + eta);
	Shape.at(2) = 0.25*(1 - ksi)*(1 - eta);
	Shape.at(3) = 0.25*(1 + ksi)*(1 - eta);

	DShape.at(0, 0) = 0.25*(1 + eta);
	DShape.at(1, 0) = 0.25*(1 + ksi);
	DShape.at(0, 1) = -0.25*(1 + eta);
	DShape.at(1, 1) = 0.25*(1 - ksi);
	DShape.at(0, 2) = -0.25*(1 - eta);
	DShape.at(1, 2) = -0.25*(1 - ksi);
	DShape.at(0, 3) = 0.25*(1 - eta);
	DShape.at(1, 3) = -0.25*(1 + ksi);

	//DShape.Print();
	//Coors->Print();
	Jacobi = DShape.Mult(Coors);
	//Jacobi.Print();
	Det = Jacobi.Determinant();

	InvJacobi = Jacobi.Inverse();
	//InvJacobi.Print();

	DShapeX = InvJacobi.Mult(DShape);
	//DShapeX.Print();
	return Jacobi;
}
FloatMatrix Quadr::ComputeBMarix(int inode)
{
	BMatrix.at(0, 0) = DShapeX.at(0, inode);
	BMatrix.at(1, 1) = DShapeX.at(1, inode);
	BMatrix.at(2, 0) = DShapeX.at(1, inode);
	BMatrix.at(2, 1) = DShapeX.at(0, inode);

	return BMatrix;
}
FloatMatrix Quadr::ComputeBMatrix()
{
	FloatMatrix BMatrixBig(3, 8);
	for (int i = 0; i < 4; i++)
	{
		BMatrixBig.at(0, i*2+0) = DShapeX.at(0, i);
		BMatrixBig.at(1, i*2+1) = DShapeX.at(1, i);
		BMatrixBig.at(2, i*2+0) = DShapeX.at(1, i);
		BMatrixBig.at(2, i*2+1) = DShapeX.at(0, i);
	}
	return BMatrixBig;
}
FloatMatrix  Quadr::ComputeConstitutiveMatrix()
{
	FloatMatrix T(3,3);
	double Young, Possion;
	Young = Mat.GetYoung();
	Possion = Mat.GetPossion();
	double D1, D2, D3;
	D1 = (1 - Possion)*Young / ((1 + Possion)*(1 - 2 * Possion));
	D2 = Possion*Young / ((1 + Possion)*(1 - 2 * Possion));
	D3 = Young / (2 * (1 + Possion));
	T.at(0, 0) = D1;
	T.at(1, 0) = D2;
	T.at(0, 1) = D2;
	T.at(1, 1) = D1;
	T.at(2, 2) = D3;
	return T;
}
FloatMatrix  Quadr::ComputeStiff()
{
	FloatArray GaussCoor(2);
	GaussPoint B;
	FloatMatrix G(2, 2);
	FloatMatrix BT(3, 2);
	FloatMatrix DT(3, 3);
	double W1,W2;
	DMatrix = ComputeConstitutiveMatrix();
	for (int inode = 0; inode < 4; inode++)
	{
		for (int jnode = 0; jnode < 4; jnode++)
		{
			for (int iksi = 0; iksi < 3; iksi++)
			{
				GaussCoor .at(0) = Gauss3[iksi];
				W1 = Weight3[iksi];
				for (int ieta = 0; ieta < 3; ieta++)
				{
					GaussCoor.at(1) = Gauss3[ieta];
					W2 = Weight3[ieta];
					B.Init(0, GaussCoor);
					ComputeJacobi(B);
					BMatrix = ComputeBMarix(inode);
					BT=BMatrix.Trans();
					//BT.Print(); 
					DT = BT.Mult(DMatrix);
					//DT.Print();
					BMatrix = ComputeBMarix(jnode);
					G = DT.Mult(BMatrix);
					W2 = W1*W2*Det;
					G = G.Mult(W2);
					//G.Print();
					Stiff.at(2 * inode, 2 * jnode) += G.at(0, 0);
					Stiff.at(2 * inode + 1, 2 * jnode) += G.at(1, 0);
					Stiff.at(2 * inode, 2 * jnode + 1) += G.at(0, 1);
					Stiff.at(2 * inode + 1, 2 * jnode + 1) += G.at(1, 1);
					
				}
			}
		}
	}
	cout << "**********************************************" << endl;
	Stiff.Print();
	return Stiff;
}
FloatArray Quadr::ComputeStress()
{
	const double a = 1 + sqrt(3) / 2, b = -1.0 / 2, c = 1 - sqrt(3) / 2;
	FloatMatrix TransMatrix(4, 4);
	TransMatrix.at(0, 0) = a;
	TransMatrix.at(0, 1) = b;
	TransMatrix.at(0, 2) = c;
	TransMatrix.at(0, 3) = b;
	TransMatrix.at(1, 0) = b;
	TransMatrix.at(1, 1) = a;
	TransMatrix.at(1, 2) = b;
	TransMatrix.at(1, 3) = c;
	TransMatrix.at(2, 0) = c;
	TransMatrix.at(2, 1) = b;
	TransMatrix.at(2, 2) = a;
	TransMatrix.at(2, 3) = b;
	TransMatrix.at(3, 0) = b;
	TransMatrix.at(3, 1) = c;
	TransMatrix.at(3, 2) = b;
	TransMatrix.at(3, 3) = a;
	Strain.Clear();
	
	GaussPoint B;
	FloatArray Coor;
	FloatArray TransTemp(4);
	Coor.SetSize(2);

	Coor.at(0) = 0;
	Coor.at(1) = 0;
	B.Init(0, Coor);
	ComputeJacobi(B);
	BMatrixBig = ComputeBMatrix();
	BMatrixBig.Print();
	Displacement.Print();
	Strain = BMatrixBig.Mult(Displacement);
	Stress = DMatrix.Mult(Strain);

	for (int inode = 0; inode < 4; inode++)
	{
		Coor.at(0) = P4ksicor[inode];
		Coor.at(1) = P4etacor[inode];
		B.Init(0, Coor);
		ComputeJacobi(B);
		BMatrixBig = ComputeBMatrix();
		GaussStrain[inode] = BMatrixBig.Mult(Displacement);
		GaussStress[inode] = DMatrix.Mult(Strain);
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			TransTemp.at(j) = GaussStrain[j].at(0);
		}
		TransTemp = TransMatrix.Mult(TransTemp);
		for (int j = 0; j < 4; j++)
		{
			NodeStrain[j].at(0) = TransTemp.at(j);
		}
		for (int j = 0; j < 4; j++)
		{
			TransTemp.at(j) = GaussStress[j].at(0);
		}
		TransTemp = TransMatrix.Mult(TransTemp);
		for (int j = 0; j < 4; j++)
		{
			NodeStress[j].at(0) = TransTemp.at(j);
		}
	}
	return NULL;
}
Quadr::~Quadr()
{
}


Quadr & Quadr::operator=(const Quadr &Q)
{
	Index = Q.Index;
	ConstitutiveMatrix = Q.ConstitutiveMatrix;
	DegreeOfFreedom = Q.DegreeOfFreedom;
	GaussPointArray = Q.GaussPointArray;
	DegreeOfFreedom = Q.DegreeOfFreedom;
	group = Q.group;
	mat = Q.mat;
	nGaussPoint = Q.nGaussPoint;
	nNodes = Q.nNodes;
	Nodes = Q.Nodes;
	Stiff = Q.Stiff;
	type = Q.type;
	return *this;
}


void Quadr::Print()
{
	cout << "Quadr Element "<<setw(5)<<this->Index;
	for (int i = 0; i < 4; i++)
	{
		cout <<setw(5)<< this->Nodes.at(i);
	}
	cout << endl;
}

Line::Line()
{
	this->type = Linear;
	this->nNodes = 2;
	this->Nodes.SetSize(2);
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
	Nodes.SetSize(nNodes);
	Nodes = L.Nodes;
	Shape.SetSize(2);
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
	Index = L.Index;
	ConstitutiveMatrix = L.ConstitutiveMatrix;
	DegreeOfFreedom = L.DegreeOfFreedom;
	GaussPointArray = L.GaussPointArray;
	DegreeOfFreedom = L.DegreeOfFreedom;
	group = L.group;
	mat = L.mat;
	nGaussPoint = L.nGaussPoint;
	nNodes = L.nNodes;
	Nodes = L.Nodes;
	Stiff = L.Stiff;
	type = L.type;
	AdjElem = L.AdjElem;
	return *this;
}

FloatArray Line::ComputeShape(GaussPoint *B)
{
	double ksi;
	ksi = B->GetCoordinate(0);
	Shape.at(0) = 0.5*(1 - ksi);
	Shape.at(1) = 0.5*(1 + ksi);
	return Shape;
}