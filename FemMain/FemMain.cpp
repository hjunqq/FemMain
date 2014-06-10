//////////////////////////////////////////////////////////////////////////  
///     COPYRIGHT NOTICE  
///     Copyright (c) 2014, 河海大学 
///     All rights reserved.  
///  
/// @file main.cpp
/// @brief  对于平面四边形单元求解
///  
///（本文件实现的功能的详述）  
///  
/// @version 0.0
/// @author 齐慧君
/// @date  2014年4月18日
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

	Fem1.ReadFiles();
	Fem2.ReadFiles();

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
	if (glb.fail() || cor.fail() || ele.fail() || grp.fail() || loa.fail() || mat.fail() || pre.fail())
	{
		cout << "Input Files Fail" << endl;
		ShowTime();
	}
	chk.open(workdir + ".chk");
	res.open(workdir + ".res");

	stream.str("");
	getline(glb, text);
	getline(glb, text);
	stream << text;
	stream >> nDim >> nNode >> nGroup >> nElem >> nMat >> nStep;

	stream.str("");
	Groups = new Group[nGroup];
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Idx, Type, nElem, Mat, Dof;
		getline(grp, text);
		getline(grp, text);
		stream << text;
		stream >> Idx >> Type >> nElem >> Mat >> Dof;
		Groups[igroup].Init(Idx, nElem, Type, Mat, Dof);
	}
	
	stream.str("");
	Nodes = new Node[nNode];
	for (int inode = 0; inode < nNode; inode++)
	{
		FloatArray Coors(nDim);
		int Idx;
		getline(cor, text);
		stream << text;
		stream >> Idx;
		Idx--;
		for (int idim = 0; idim < nDim; idim++)
		{
			stream >> Coors.at(idim);
		}
		Nodes[inode].Init(Idx, &Coors);
	}
	for (int inode = 0; inode < nNode; inode++)
	{
		Nodes[inode].Print();
	}
}