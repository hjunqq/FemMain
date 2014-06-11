//////////////////////////////////////////////////////////////////////////  
///     COPYRIGHT NOTICE  
///     Copyright (c) 2014, �Ӻ���ѧ 
///     All rights reserved.  
///  
/// @file main.cpp
/// @brief  ����ƽ���ı��ε�Ԫ���
///  
///�����ļ�ʵ�ֵĹ��ܵ�������  
///  
/// @version 0.0
/// @author ��۾�
/// @date  2014��4��18��
///  
///  
/// 
//////////////////////////////////////////////////////////////////////////

#include "FemMain.h"
int main(int argc,char**argv)
{
	FemMain Fem1;
	FemMain Fem2;
	string text;
	stringstream stream;
	ifstream inp("1.inp");
	getline(inp, text);

	getline(inp, text);
	Fem1.WordDir()=text;
	getline(inp, text);
	Fem2.WordDir()=text;
	
	inp.close();
	Fem1.ShowTime();

	Fem1.ReadFiles();
	Fem2.ReadFiles();

	Fem1.GIDOutMesh();

	Fem1.ShowTime();
	return 0;
}
void FemMain::ShowTime()
{
	time(&tNow);
	localtime_s(&tmLocal, &tNow);
	cout << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon + 1 << "-" <<
		tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0') << setw(2) << tmLocal.tm_min << ":"
		<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
	chk << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon << "-" <<
		tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0') << setw(2) << tmLocal.tm_min << ":"
		<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
}
string & FemMain::WordDir()
{
	return workdir;
}
void FemMain::ReadFiles()
{
	stringstream stream;
	glb.open(workdir + ".glb");
	cor.open(workdir + ".cor");
	ele.open(workdir + ".ele");
	grp.open(workdir + ".grp");
	loa.open(workdir + ".loa");
	mat.open(workdir + ".mat");
	pre.open(workdir + ".pre");

	chk.open(workdir + ".chk");

	if (glb.fail() || cor.fail() || ele.fail() || grp.fail() || loa.fail() || mat.fail() || pre.fail())
	{
		cout << "Input Files Fail" << endl;
		ShowTime();
	}
	chk.open(workdir + ".chk");
	res.open(workdir + ".res");

	stream.clear();
	stream.str("");
	glb.getline(str, MAXCHAR);
	glb.getline(str, MAXCHAR);
	stream << str;
	stream >> nDim >> nNode >> nGroup >> nElem >> nMat >> nStep;


	Groups = new Group[nGroup];
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Idx, Type, nElem, Mat, Dof;
		grp.getline(str, MAXCHAR);

		stream.clear();
		stream.str("");
		grp.getline(str, MAXCHAR);
		stream << str;
		stream >> Idx >> Type >> nElem >> Mat >> Dof;
		Groups[igroup].Init(Idx, nElem, Type, Mat, Dof);
	}

	Nodes = new Node[nNode];
	FloatArray *Coors;
	Coors = new FloatArray(nDim);
	for (int inode = 0; inode < nNode; inode++)
	{
		int Idx;
		cor.getline(str, MAXCHAR);
		stream.clear();
		stream.str("");
		stream << str;
		stream >> Idx;
		Idx--;
		for (int idim = 0; idim < nDim; idim++)
		{
			stream >> Coors->at(idim);
		}
		Nodes[inode].Init(Idx, Coors);
	}
	for (int inode = 0; inode < nNode; inode++)
	{
		Nodes[inode].Print();
	}

	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int type = Groups[igroup].GetType();
		int nElem = Groups[igroup].GetnElements();
		int Mat = Groups[igroup].GetMaterial();
		if (type == Quadrilateral)
		{
			int Idx;
			IntArray Nodes(4);
			Quadr *Quadrs = new Quadr[nElem]();
			for (int ielem = 0; ielem < nElem; ielem++)
			{
				stream.clear();
				stream.str("");
				ele.getline(str, MAXCHAR);
				stream << str;
				stream >> Idx;
				for (int i = 0; i < 4; i++)
				{
					stream >> Nodes.at(i);
					Nodes.at(i)--;
				}
				Quadrs[ielem].Init(4, Mat, Idx, igroup, &Nodes);
			}
			for (int ielem = 0; ielem < nElem; ielem++)
			{
				Quadrs[ielem].Print();
			}
		}
	}


	Mats = new Material[nMat]();
	for (int imat = 0; imat < nMat; imat++)
	{
		int Idx, Type;
		Idx = imat;
		mat.getline(str, MAXCHAR);
		if (strcmp(str, "Elastic") == 0)
		{
			Type = Elastic;
			double Young, Possion, Density;
			stream.clear();
			stream.str("");
			mat.getline(str, MAXCHAR);
			stream << str;
			stream >> Young >> Possion >> Density;
			Mats[imat].Init(Idx, Elastic, Young, Possion, Density);
		}
	}

	loa.getline(str, MAXCHAR);
	stream.clear();
	stream.str("");
	loa.getline(str, MAXCHAR);
	stream << str;
	stream >> nTotalLoad >> nFace >> nConcentrate >> nVolumn;
	Faces = new Face[nFace]();
	Cons = new Concentrate[nConcentrate]();
	Vols = new Volumn[nVolumn]();
	for (int iface = 0; iface < nFace; iface++)
	{
		int nEle = 0, nNode = 0, Dir = 0, AdjElem = 0;
		double StartC, EndC, StartV, EndV;
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);

		stream.str("");
		stream.clear();
		stream << str;
		stream >> nEle >> nNode >> Dir >> StartC >> EndC >> StartV >> EndV;
		Faces[iface].Init(iface, nEle, nNode, StartC, StartV, EndC, EndV, Dir);
		loa.getline(str, MAXCHAR);
		if (nNode == 2)
		{
			Line *Lines;
			Lines = new Line[nEle]();
			IntArray Nodes(2);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				loa.getline(str, MAXCHAR);
				stream.clear();
				stream.str("");
				stream << str;
				stream >> Nodes.at(0) >> Nodes.at(1) >> AdjElem;
				Nodes.at(0)--;
				Nodes.at(1)--;
				AdjElem--;
				Lines[ielem].Init(2, NULL, ielem, NULL, &Nodes);
				Lines[ielem].AtAdjElem() = AdjElem;
			}
			Faces[iface].Set(Lines);
			delete[] Lines;
		}
	}
	for (int ivol = 0; ivol < nVolumn; ivol++)
	{
		int Dir,nGroup;
		double Acc;
		IntArray *group;
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> Dir >> Acc >> nGroup;
		group=new IntArray(nGroup);
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		for (int igroup = 0; igroup < nGroup; igroup++)
		{
			stream >> group->at(igroup);
			group->at(igroup)--;
		}
		Vols[ivol].Init(ivol, nGroup, group, Acc, Dir);
		delete[] group;
	}
	for (int icon = 0; icon < nConcentrate; icon++)
	{
		int nNode, Dir;
		FloatArray *Value;
		IntArray *Node;
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> nNode >> Dir;
		Value = new FloatArray(nNode);
		Node = new IntArray(nNode);
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Value->at(inode);
		}
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		stream.str ("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Node->at(inode);
			Node->at(inode)--;
		}
		Cons[icon].Init(icon, nNode, Node, Value, Dir);
	}

	
	pre.getline(str, MAXCHAR);
	pre.getline(str, MAXCHAR);
	stream.str("");
	stream.clear();
	stream << str;
	stream >>nTotalBoundary>> nFix >> nDisp >> nInter;
	Fix = new Fixed[nFix];
	Disp = new Displace[nDisp];
	Inters = new Interact[nInter];
	for (int ifix = 0; ifix < nFix; ifix++)
	{
		int Dir, nNode;
		IntArray *Nodes;
		pre.getline(str, MAXCHAR);
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> Dir >> nNode;
		Nodes = new IntArray(nNode);
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Nodes->at(inode);
			Nodes->at(inode)--;
		}
		Fix[ifix].Init(ifix, nNode, Nodes, Dir);
	}
	for (int idisp = 0; idisp < nDisp; idisp++)
	{
		int Dir, nNode;
		IntArray *Node;
		FloatArray *Value;
		pre.getline(str, MAXCHAR);
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> Dir >> nNode;
		Node = new IntArray(nNode);
		Value = new FloatArray(nNode);
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Value->at(inode);
		}
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Node->at(inode);
			Node->at(inode)--;
		}
		Disp[idisp].Init(nNode, Dir, Node, Value);
	}
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		int nNode, AdjDomain;
		IntArray *Local, *Remote;
		pre.getline(str, MAXCHAR);
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> nNode >> AdjDomain;
		Local = new IntArray(nNode);
		Remote = new IntArray(nNode);
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Local->at(inode);
			Local->at(inode)--;
		}
		pre.getline(str, MAXCHAR);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Remote->at(inode);
			Remote->at(inode)--;
		}
		Inters[iinter].Init(nNode, AdjDomain, Local, Remote);
	}

	glb.close();
	cor.close();
	ele.close();
	grp.close();
	loa.close();
	mat.close();
	pre.close();
}

void FemMain::GIDOutMesh()
{
	string GidMeshFile;
	FloatArray *Coor;
	int *ENode;
	GidMeshFile = workdir + ".flavia.msh";
	GiD_OpenPostMeshFile(GidMeshFile.c_str(), GiD_PostAscii);
	GiD_BeginCoordinates();
	Coor = new FloatArray(nDim);
	for (int inode = 0; inode < nNode; inode++)
	{
		*Coor = Nodes->GetCoordinate();
		if (nDim == 2)
		{
			GiD_WriteCoordinates2D(inode + 1, Coor->at(0), Coor->at(1));
		}

	}
	GiD_EndCoordinates();
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int type,nEle,nNode;
		type = Groups->GetType();
		nEle = Groups->GetnElements();
		GiD_BeginMeshGroup(GroupName[type].c_str);
		for (int iele = 0; iele < nEle; iele++)
		{
		}
	}
	GiD_BeginElements();

	GiD_EndElements();
	GiD_EndMesh();
	GiD_ClosePostMeshFile();
}