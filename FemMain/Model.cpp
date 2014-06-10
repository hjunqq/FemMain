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


 // 单元初始化
void Element::Init(int nNodes, int Material, int Index, int Group, IntArray Node)
{
	this->nNodes = nNodes;
	this->Material = Material;
	this->Index = Index;
	this->group = Group;
	this->Nodes = Node;
	this->DegreeOfFreedom = NULL;
	this->Stiff = NULL;
}


// 计算刚度矩阵
FloatMatrix * Element::BuildStiff()
{
	return NULL;
}


// 组装本构矩阵
FloatMatrix * Element::BuildConstitutiveMatrix()
{
	return NULL;
}


// 计算高斯点应变
FloatArray * Element::ComputeStrain(GaussPoint * B)
{
	return NULL;
}


// 计算B矩阵
FloatMatrix Element::ComputeBMarix(GaussPoint * B)
{
	return FloatMatrix();
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
IntArray Element::GetNodeArray()
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
	return Material;
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
	this->Coordinates = Coordinate;
}


// 打印节点信息
void Node::Print()
{
	cout << "Node " << Index << endl;
	Coordinates->Print();
}


// 获得节点坐标
double Node::GetCoordinate(int i)
{
	return Coordinates->at(i);
}


// 获得节点坐标
FloatArray Node::GetCoordinate()
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
	return Stress;
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
void Group::Init(int Index, int nElements, int type)
{
	this->Index = Index;
	this->nElements = nElements;
	this->type = type;
	Elements = new IntArray(nElements);
}


// 填充单元
void Group::FillElement(IntArray * ElementList)
{
	*Elements = Elements->Copy(ElementList);
}


// 获取第i个单元号
int Group::GetElement(int i)
{
	return Elements->at(i);
}




// 获取单元个数
int Group::GetnElements()
{
	return nElements;
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
}


Quadr::~Quadr()
{
}


void Element::Print()
{
}


void Quadr::Print()
{
	cout << "I am a Quadr Element" << endl;
}

Line::Line()
{
	this->type = Linear;
}
Line::~Line()
{
}


void Line::Print()
{
	cout << "I am a Line Element" << endl;
}
