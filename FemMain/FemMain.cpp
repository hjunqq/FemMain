//////////////////////////////////////////////////////////////////////////  
///     COPYRIGHT NOTICE  
///     Copyright (c) 2014, 河海大学 
///     All rights reserved.  
///  
/// @file main.cpp
/// @brief  对于平面四边形单元静力求解，模态计算
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
	//_CrtSetBreakAlloc(65741);
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	//_CrtMemState s1;
	//_CrtMemState s2;
	//_CrtMemCheckpoint(&s1);
	FemMain Fem;
	string text1,text2;
	stringstream stream;
	ifstream inp("1.inp");
	int nProb;

	//bool *Converge;
	//double Error;

	int NProces, MyID;
	MPI::Init(argc, argv);
	NProces = MPI::COMM_WORLD.Get_size();
	MyID = MPI::COMM_WORLD.Get_rank();
	Fem.GetSize(NProces);
	Fem.GetID(MyID);
	cout << setw(10) << NProces << setw(10) << MyID << endl;
	getline(inp, text1);
	getline(inp, text1);
	stream << text1;
	stream >> nProb;
	for (int iProb = 0; iProb < nProb; iProb++)
	{
		getline(inp, text1);
		if (MyID == iProb)
		{
			Fem.WorkDir() = text1;
		}
	}
	//getline(inp, text1);
	//getline(inp, text2);
	//if (MyID == 0)
	//{
	//	Fem.WorkDir() = text1;
	//}
	//else if (MyID == 1)
	//{
	//	Fem.WorkDir() = text2;
	//}
	//
	inp.close();
	MPI::COMM_WORLD.Barrier();
	GiD_PostInit();

	Fem.ShowTime();

	Fem.ReadFiles();

	Fem.OpenGidFile(); 


	Fem.GIDOutMesh();

	Fem.ShowTime();
	MPI::COMM_WORLD.Barrier();

	Fem.GloableSolve();


	Fem.CloseGidFile();

	GiD_PostDone();
	MPI::COMM_WORLD.Barrier();
	MPI::Finalize();
	//_CrtMemCheckpoint(&s2);
	//_CrtMemState s3;
	//if (_CrtMemDifference(&s3, &s1, &s2))
	//{
	//	_CrtMemDumpStatistics(&s3);
	//}
	//_CrtDumpMemoryLeaks();
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
string & FemMain::WorkDir()
{
	return workdir;
}
int FemMain::GetStep()
{
	return nStep;
}
int FemMain::GetMaxIter()
{
	return MaxIter;
}
void FemMain::ReadFiles()
{
	stringstream stream;
	cout << workdir << setw(10) << Id << endl;
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
		cout << "Input Files Fail" <<setw(10)<<Id<< endl;
		ShowTime();
	}

	

	stream.clear();
	stream.str("");
	glb.getline(str, MACLENGTH);
	glb.getline(str, MACLENGTH);
	stream << str;
	stream >> nDim >> nNode >> nGroup >> nElem >> nMat >> nStep >> nDof >> MaxIter >> Tolerance;
	glb.getline(str, MACLENGTH);
	glb.getline(str, MACLENGTH);
	stream.clear();
	stream.str("");
	stream << str;
	stream >> SolveMethod >> SolveType>>ProblemType>>dT;
	cout << Tolerance << endl;

	Groups = new Group[nGroup];
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Idx, Type, nElem, Mat, Dof;
		double alpha, beta;
		grp.getline(str, MACLENGTH);

		stream.clear();
		stream.str("");
		grp.getline(str, MACLENGTH);
		stream << str;
		stream >> Idx >> Type >> nElem >> Mat >> Dof;
		Mat--;
		stream.clear();
		stream.str("");
		grp.getline(str, MACLENGTH);
		grp.getline(str, MACLENGTH);
		stream << str;
		stream >> alpha >> beta;
		Groups[igroup].Init(Idx, nElem, Type, Mat, Dof);
		Groups[igroup].FillDampPara(alpha, beta);
	}

	Nodes = new Node[nNode];
	FloatArray Coors;
	Coors.SetSize(nDim);
	for (int inode = 0; inode < nNode; inode++)
	{
		int Idx;
		cor.getline(str, MACLENGTH);
		stream.clear();
		stream.str("");
		stream << str;
		stream >> Idx;
		Idx--;
		for (int idim = 0; idim < nDim; idim++)
		{
			stream >> Coors.at(idim);
		}
		Nodes[inode].Init(Idx, Coors);
	}
	//for (int inode = 0; inode < nNode; inode++)
	//{
	//	Nodes[inode].Print();
	//}

	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int type = Groups[igroup].GetType();
		int nElem = Groups[igroup].GetnElements();
		int Mat = Groups[igroup].GetMaterial();
		int Dof = Groups[igroup].GetDof();
		if (type == Quadrilateral)
		{
			int Idx;
			IntArray Nodes(4);
			Quadr *Quadrs = new Quadr[nElem];
			for (int ielem = 0; ielem < nElem; ielem++)
			{
				stream.clear();
				stream.str("");
				ele.getline(str, MACLENGTH);
				stream << str;
				stream >> Idx;
				for (int i = 0; i < 4; i++)
				{
					stream >> Nodes.at(i);
					Nodes.at(i)--;
				}
				Quadrs[ielem].Init(4, Mat, Idx, igroup,Dof, &Nodes);
			}
			//for (int ielem = 0; ielem < nElem; ielem++)
			//{
			//	Quadrs[ielem].Print();
			//}
			Groups[igroup].FillElement(Quadrs);
			delete[] Quadrs;
			Quadrs = NULL;
		}
	}


	Mats = new Material[nMat]();
	for (int imat = 0; imat < nMat; imat++)
	{
		int Idx, Type;
		Idx = imat;
		mat.getline(str, MACLENGTH);
		if (strcmp(str, "Elastic") == 0)
		{
			Type = Elastic;
			double Young, Possion, Density;
			stream.clear();
			stream.str("");
			mat.getline(str, MACLENGTH);
			stream << str;
			stream >> Young >> Possion >> Density;
			Mats[imat].Init(Idx, Elastic, Young, Possion, Density);
		}
	}

	loa.getline(str, MACLENGTH);
	stream.clear();
	stream.str("");
	loa.getline(str, MACLENGTH);
	stream << str;
	stream >> nTotalLoad >> nFace >> nVolumn >> nConcentrate;
	Faces = new Face[nFace]();
	Cons = new Concentrate[nConcentrate]();
	Vols = new Volumn[nVolumn]();
	for (int iface = 0; iface < nFace; iface++)
	{
		int nEle = 0, nNode = 0, Dir = 0,AdjElem=0;
		double StartC, EndC, StartV, EndV;
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);

		stream.str("");
		stream.clear();
		stream << str;
		stream >> nEle >> nNode >> Dir >> StartC >> EndC >> StartV >> EndV;
		Dir--;
		Faces[iface].Init(iface, nEle, nNode, StartC, StartV, EndC, EndV, Dir);
		loa.getline(str, MACLENGTH);
		if (nNode == 2)
		{
			Line *Lines;
			Lines = new Line[nEle]();
			IntArray Nodes(2);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				loa.getline(str, MACLENGTH);
				stream.clear();
				stream.str("");
				stream << str;
				stream >> Nodes.at(0) >> Nodes.at(1) >> AdjElem;
				Nodes.at(0)--;
				Nodes.at(1)--;
				AdjElem--;
				Lines[ielem].Init(2, NULL, ielem, NULL,NULL, &Nodes);
				Lines[ielem].AtAdjElem() = AdjElem;
			}
			Faces[iface].Set(Lines);
			delete[] Lines;
			Lines = NULL;
		}
	}
	for (int ivol = 0; ivol < nVolumn; ivol++)
	{
		int Dir,nGroup;
		double Acc;
		IntArray group;
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> Dir >> Acc >> nGroup;
		Dir--;
		group.SetSize(nGroup);
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		for (int igroup = 0; igroup < nGroup; igroup++)
		{
			stream >> group.at(igroup);
			group.at(igroup)--;
		}
		Vols[ivol].Init(ivol, nGroup, group, Acc, Dir);
	}
	for (int icon = 0; icon < nConcentrate; icon++)
	{
		int nNode, Dir;
		FloatArray Value;
		IntArray Node;
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> nNode >> Dir;
		Dir--;
		Value.SetSize(nNode);
		Node.SetSize(nNode);
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Value.at(inode);
		}
		loa.getline(str, MACLENGTH);
		loa.getline(str, MACLENGTH);
		stream.str ("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Node.at(inode);
			Node.at(inode)--;
		}
		Cons[icon].Init(icon, nNode,Node,Value, Dir);
	}

	
	pre.getline(str, MACLENGTH);
	pre.getline(str, MACLENGTH);
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
		IntArray Nodes;
		pre.getline(str, MACLENGTH);
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> Dir >> nNode;
		Dir--;
		Nodes.SetSize(nNode);
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Nodes.at(inode);
			Nodes.at(inode)--;
		}
		Fix[ifix].Init(ifix, nNode,Nodes, Dir);
	}
	for (int idisp = 0; idisp < nDisp; idisp++)
	{
		int Dir, nNode;
		IntArray Node;
		FloatArray Value;
		pre.getline(str, MACLENGTH);
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> Dir >> nNode;
		Dir--;
		Node.SetSize(nNode);
		Value.SetSize(nNode);
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Value.at(inode);
		}
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Node.at(inode);
			Node.at(inode)--;
		}
		Disp[idisp].Init(nNode, Dir, Node, Value);
	}
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		int nNode, AdjDomain;
		IntArray Local, Remote;
		pre.getline(str, MACLENGTH);
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		stream >> nNode >> AdjDomain;
		AdjDomain--;
		Local.SetSize(nNode);
		Remote.SetSize(nNode);
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Local.at(inode);
			Local.at(inode)--;
		}
		pre.getline(str, MACLENGTH);
		stream.str("");
		stream.clear();
		stream << str;
		for (int inode = 0; inode < nNode; inode++)
		{
			stream >> Remote.at(inode);
			Remote.at(inode)--;
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
void FemMain :: OpenGidFile()
{
	string GidMeshFile, GidResultFile;
	
	GidMeshFile = workdir + ".flavia.msh";
	GidResultFile = workdir + ".flavia.res";
	FMesh = GiD_fOpenPostMeshFile(GidMeshFile.c_str(), GiD_PostAscii);
	FRes=GiD_fOpenPostResultFile(GidResultFile.c_str(), GiD_PostAscii);
}
void FemMain::GIDOutMesh()
{
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int type,nEle;
		int *nid;
		type = Groups[igroup].GetType();
		nEle = Groups[igroup].GetnElements();
		GiD_fBeginMesh(FMesh, GroupName[type - 1].c_str(), GiD_Dimension(nDim), GiD_ElementType(type), type);
		if (igroup == 0)
		{
			FloatArray Coor;
			GiD_fBeginCoordinates(FMesh);
			Coor.SetSize(nDim);
			for (int inode = 0; inode < nNode; inode++)
			{
				Coor = Nodes[inode].GetCoordinate();
				if (nDim == 2)
				{
					GiD_fWriteCoordinates2D(FMesh,inode + 1, Coor.at(0), Coor.at(1));
				}
				else
				{
					GiD_fWriteCoordinates(FMesh,inode + 1, Coor.at(0), Coor.at(1), Coor.at(2));
				}
			}
			GiD_fEndCoordinates(FMesh);
		}
		GiD_fBeginElements(FMesh);
		if (type == 4)
		{
			Quadr *Elems=NULL;
			IntArray Nodes;
			nid = new int[type];
			Elems =  Groups[igroup].GetElement(*Elems);
			for (int iele = 0; iele < nEle; iele++)
			{
				Nodes = Elems[iele].GetNodeArray();
				for (int inode = 0; inode < type; inode++)
				{
					nid[inode] = Nodes.at(inode)+1;
				}
				GiD_fWriteElement(FMesh,Elems[iele].GetIndex(), nid);
			}
		}
		
		GiD_fEndElements(FMesh);
		GiD_fEndMesh(FMesh);
	}
	GiD_fClosePostMeshFile(FMesh);
}

void FemMain::ComputeDOF()
{
	DegreeOfFreedom.SetSize(nDof*nNode);
	DegreeOfFreedom.Set(1);
	for (int ifix = 0; ifix < nFix; ifix++)
	{
		int iDir = Fix[ifix].GetDir();
		int nFixNode = Fix[ifix].GetnNode();
		int FixNode;
		for (int inode = 0; inode < nFixNode; inode++)
		{
			FixNode = Fix[ifix].at(inode);
			DegreeOfFreedom.at(FixNode*nDof + iDir) = 0;
		}
	}
	for (int idisp = 0; idisp < nDisp; idisp++)
	{
		int iDir = Disp[idisp].GetDir();
		int nDispNode = Disp[idisp].GetnNode();
		int DispNode;
		for (int inode = 0; inode < nDispNode; inode++)
		{
			DispNode = Disp[idisp].at(inode);
			DegreeOfFreedom.at(DispNode*nDof + iDir) = 0;
		}
	}
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		int nInterNode = Inters[iinter].GetnNode();
		IntArray Local(Inters[iinter].GetLocal());
		int InterNode;
		for (int inode = 0; inode < nInterNode; inode++)
		{
			InterNode = Local.at(inode);
			for (int iDof = 0; iDof < nDof; iDof++)
			{
				DegreeOfFreedom.at(InterNode*nDof + iDof) = 0;
			}
		}
	}
	TotalDOF = 0;
	for (int idof = 0; idof < nDof*nNode; idof++)
	{
		if (DegreeOfFreedom.at(idof))
		{
			TotalDOF++;
			DegreeOfFreedom.at(idof) = TotalDOF;
		}
	}
	//cout << "DegreeOfFreedom      " ;
	//DegreeOfFreedom.Print();
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems=NULL;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem].FillDof(DegreeOfFreedom);
			}
		}
	}
}

void FemMain::InitSolve()
{
	Stiff.SetSize(TotalDOF, TotalDOF);
	Mass.SetSize(TotalDOF, TotalDOF);
	ResultZero .SetSize(TotalDOF);
	LResultZero.SetSize(TotalDOF);
	ResultFirst.SetSize(TotalDOF);
	LResultFirst.SetSize(TotalDOF);
	ResultSecond.SetSize(TotalDOF);
	LResultSecond.SetSize(TotalDOF);
	ExternalForce .SetSize(TotalDOF);
	InitialStain .SetSize(TotalDOF);
	InitialDispLoad .SetSize(TotalDOF);
	IniDisplacement.SetSize(nNode*nDof);
	InterDisplace.SetSize(nNode*nDof);
	EffictiveLoad.SetSize(TotalDOF);
	TotalLoad.SetSize(TotalDOF);
	InteractLoad.SetSize(TotalDOF);

	InteractValue = new FloatArray[nInter];
	InteractValueOld = new FloatArray[nInter];
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		InteractNode = GetInteractNode(iinter);
		int size = InteractNode.GetSize();
		InteractValue[iinter].SetSize(size*nDof);
		InteractValueOld[iinter].SetSize(size*nDof);
	}

	Eigenvalues.SetSize(MaxIter);
	Converge = new bool[2];
	Converge[0] = false;
	Converge[1] = false;
	int NonZero;
	IntArray *RowIdx, *ColIdx,*CNonZero;
	RowIdx = new IntArray(TotalDOF+1);
	CNonZero = new IntArray(TotalDOF);

	IntMatrix *DofRel;
	DofRel = new IntMatrix(TotalDOF,TotalDOF);
	
	IntArray *EDof;
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		int Dof = Groups[igroup].GetDof();
		if (Type == 4)
		{
			Quadr *Elems=NULL;
			EDof = new IntArray(Type*Dof);
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				*EDof=Elems[ielem].GetDof();
				for (int iDof = 0; iDof < Type*Dof; iDof++)
				{
					if (EDof->at(iDof) !=0 )
					{
						for (int jDof = 0; jDof < Type*Dof; jDof++)
						{
							if (EDof->at(jDof) != 0)
							{
								DofRel->at(EDof->at(iDof) - 1, EDof->at(jDof) - 1) = 1;
							}
						}
					}
				}
			}
		}
	}
	NonZero = 0;
	for (int iDof = 0; iDof < TotalDOF; iDof++)
	{
		for (int jDof = 0; jDof < TotalDOF; jDof++)
		{
			CNonZero->at(iDof) += DofRel->at(iDof, jDof);
		}
	}
	RowIdx->at(0) = 0;
	for (int iDof = 0; iDof < TotalDOF; iDof++)
	{
		RowIdx->at(iDof + 1) = RowIdx->at(iDof) + CNonZero->at(iDof);
	}
	NonZero = RowIdx->at(TotalDOF);
	ColIdx = new IntArray(NonZero);
	int iRow,iCol;
	for (int iDof = 0; iDof < TotalDOF; iDof++)
	{
		iCol = 0;
		iRow = RowIdx->at(iDof);
		for (int jDof = 0; jDof < TotalDOF; jDof++)
		{
			if (DofRel->at(iDof, jDof) !=0)
			{
				iCol++;
				ColIdx->at(iRow+iCol-1)=jDof;
			}
		}
	}
	//ColIdx->Print();
	
}

void FemMain::ComputeElementStiff()
{
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems=NULL;
			IntArray ENode(Type);
			FloatMatrix Coor;
			Material Mat;
			Coor.SetSize(Type,nDim);
			Mat =Mats[Groups[igroup].GetMaterial()];
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ENode = Elems[ielem].GetNodeArray();
				Elems[ielem].SetMaterial(Mat);
				for (int inode = 0; inode < Type; inode++)
				{
					for (int idim = 0; idim < nDim; idim++)
					{
						Coor.at(inode, idim) = Nodes[ENode.at(inode)].GetCoordinate().at(idim);
					}
				}
				Elems[ielem].SetCoor(Coor);
				Elems[ielem].ComputeStiff();
				Elems[ielem].ComputeMassMatrix();
			}
		}
	}
}

void FemMain::AssembleStiff()
{
	int iedof, jedof;
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems=NULL;
			FloatMatrix EStiff(Type*nDof, Type*nDof),EMass(Type*nDof,Type*nDof);
			Elems = Groups[igroup].GetElement(*Elems);
			IntArray Dof(Type*nDof);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Dof = Elems[ielem].GetDof();
				EStiff = Elems[ielem].GetStiff();
				EMass = Elems[ielem].GetMass();
				for (int idof = 0; idof < Type*nDof; idof++)
				{
					iedof = Dof.at(idof);
					for (int jdof = 0; jdof < Type*nDof; jdof++)
					{
						jedof = Dof.at(jdof);
						if (iedof != 0 && jedof != 0)
						{
							Stiff.at(iedof - 1, jedof - 1) += EStiff.at(idof, jdof);
							Mass.at(iedof - 1, jedof - 1) += EMass.at(idof, jdof);
						}
					}
				}
			}
			//delete[]Elems;
		}
	}
	//cout << "StiffMatrix=:";
	//Stiff.Print();
	//cout << "MassMatrix=:";
	//Mass.Print();

	if (DynamicStatic == ProblemType)
	{
		double Alpha, Beta;
		Alpha = Groups[0].GetDampAlpha();
		Beta = Groups[0].GetDampBeta();
		Damp.SetSize(TotalDOF, TotalDOF);
		for (int i = 0; i < TotalDOF; i++)
		{
			for (int j = 0; j < TotalDOF; j++)
			{
				Damp.at(i, j) = Mass.at(i, j)*Alpha + Stiff.at(i, j)*Beta;
			}
		}
		//Damp = Mass.Mult(Alpha) + Stiff.Mult(Beta);
		//cout << "DampMatrix=:";
		//Damp.Print();
	}
}

void FemMain::ApplyLoad()
{
	for (int iface = 0; iface < nFace; iface++)
	{
		int nEle = 0, nNode = 0, Dir = 0;
		double StartC, EndC, StartV, EndV;
		
		nEle = Faces[iface].GetnElem();
		nNode = Faces[iface].GetnNode();
		Dir = Faces[iface].GetDir();
		StartC = Faces[iface].GetStartC();
		EndC = Faces[iface].GetEndC();
		StartV = Faces[iface].GetStartV();
		EndV = Faces[iface].GetEndV();
		if (nNode == 2)
		{
			Line *Lines=NULL;
			int Type = 2;
			FloatMatrix Coor;
			IntArray ENode;
			ENode.SetSize(2);
			double Length;
			FloatArray PVal, Normal, Shape, GaussCoor, TELoad, ELoad;
			PVal.SetSize(2);
			Normal.SetSize(3);
			Shape.SetSize(2);
			ELoad.SetSize(2);
			TELoad.SetSize(2);
			GaussCoor.SetSize(1);
			GaussPoint B;
			double Det,PLoad; 
			Lines = Faces[iface].GetLines();

			Coor.SetSize(Type, nDim);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ELoad.Clear();
				ENode = Lines[ielem].GetNodeArray();
				
				for (int inode = 0; inode < Type; inode++)
				{
					for (int idim = 0; idim < nDim; idim++)
					{
						Coor.at(inode, idim) = Nodes[ENode.at(inode)].GetCoordinate().at(idim);
					}
					PVal.at(inode) = (Coor.at(inode, Dir) - StartC) / (EndC - StartC)*(EndV - StartV) + StartV;
				}
				Length = sqrt(pow(Coor.at(0, 0) - Coor.at(1, 0), 2) +
					pow(Coor.at(0, 1) - Coor.at(1, 1), 2));
				Normal.at(0) = (Coor.at(0, 0) - Coor.at(1, 0)) / Length;
				Normal.at(1) = (Coor.at(0, 1) - Coor.at(1, 1)) / Length;
				FloatArray Descartes(3);
				Descartes.at(2) = 1;
				Normal = Normal.Cross(Descartes);
				double ksi, weight;
				for (int iksi = 0; iksi < 2; iksi++)
				{
					ksi = Gauss2[iksi];
					weight = Weight2[iksi];
					Shape.at(0) = 0.5*(1 - ksi);
					Shape.at(1) = 0.5*(1 + ksi);
					Det = 0.5*Length;
					PLoad = PVal.Dot(&Shape);
					PLoad = PLoad*Det*weight;
					TELoad = Shape.Times(PLoad);
					ELoad.at(0) += TELoad.at(0);
					ELoad.at(1) += TELoad.at(1);
				}
				for (int inode = 0; inode < Type; inode++)
				{
					int NodeIdx = ENode.at(inode);
					for (int idof = 0; idof < nDof; idof++)
					{
						if (DegreeOfFreedom.at(NodeIdx * 2 + idof) != 0)
						{
							ExternalForce.at(DegreeOfFreedom.at(NodeIdx * 2 + idof) - 1) -= ELoad.at(idof)*Normal.at(idof);
						}
					}
				}
			}
		}
	}
	for (int iCon = 0; iCon < nConcentrate; iCon++)
	{
		int nNode,Dir;
		nNode = Cons[iCon].GetnNode();
		FloatArray Values(nNode);
		IntArray Nodes(nNode);
		Values = Cons[iCon].GetValues();
		Nodes = Cons[iCon].GetNodes();
		Dir = Cons[iCon].GetDir();
		for (int inode = 0; inode < nNode; inode++)
		{
			if (DegreeOfFreedom.at(Nodes.at(inode)*nDof + Dir) != 0)
			{
				ExternalForce.at(DegreeOfFreedom.at(Nodes.at(inode)*nDof + Dir)-1) += Values.at(inode);
			}
		}
	}
	for (int iVol = 0; iVol < nVolumn; iVol++)
	{
		int nBodyGroup, Dir,Type,nEle;
		int GroupIdx, iMat;
		double Density;
		double Value;
		nBodyGroup = Vols[iVol].GetnGroup();
		Dir = Vols[iVol].GetDir();
		IntArray VolGroups(nBodyGroup);
		VolGroups= Vols[iVol].GetGroup();
		Value = Vols[iVol].GetValue();
		for (int igroup = 0; igroup < nBodyGroup; igroup++)
		{
			
			GroupIdx = VolGroups.at(igroup);
			nEle = Groups[GroupIdx].GetnElements();
			Type = Groups[GroupIdx].GetType();
			iMat = Groups[GroupIdx].GetMaterial();
			Density = Mats[iMat].GetDensity();
			if (Type == 4)
			{
				Quadr *Elems=NULL;
				Elems = Groups[GroupIdx].GetElement(*Elems);
				FloatArray *Eload = new FloatArray(2);
				FloatArray GaussT(2);
				FloatArray Shape(Type);
				FloatArray TLoad(Type);
				IntArray Nodes(Type);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					TLoad.Clear();
					Nodes = Elems[ielem].GetNodeArray();
					for (int iksi = 0; iksi < 3; iksi++)
					{
						double ksi = Gauss3[iksi];
						double W1 = Weight3[iksi];
						GaussT.at(0) = ksi;
						for (int ieta = 0; ieta < 3; ieta++)
						{
							double eta = Gauss3[ieta];
							double W2 = Weight3[ieta];
							GaussT.at(2) = eta;
							GaussPoint B;
							double Det;
							B.Init(0, GaussT);
							Elems[ielem].ComputeJacobi(B);
							Shape = Elems[ielem].GetShape();
							Det = Elems[ielem].GetDet();
							for (int inode = 0; inode < Type; inode++)
							{
								TLoad.at(inode) += Shape.at(inode)*Density*Value*Det*W1*W2;
							}
						}
					}
					for (int inode = 0; inode < Type; inode++)
					{
						if (DegreeOfFreedom.at(Nodes.at(inode)*nDof + Dir) != 0)
						{
							ExternalForce.at(DegreeOfFreedom.at(Nodes.at(inode)*nDof + Dir) - 1) += TLoad.at(inode);
						}
					}
				}
			}
		}
		
	}
	for (int iDisp = 0; iDisp < nDisp; iDisp++)
	{
		FloatArray NodeDisplacement;
		FloatMatrix A;
		FloatMatrix LStiff,EStiff;
		int iedof, jedof;
		IntArray EDof;
		int NodeIndex, iDir, nDispDof, nDispNode;
		IntArray DispDof;
		IntArray DispNode;

		A.SetSize(TotalDOF, TotalDOF);
		
		DispDof = DegreeOfFreedom;
		nDispDof = TotalDOF;
		iDir = Disp[iDisp].GetDir();
		nDispNode = Disp[iDisp].GetnNode();
		
		DispNode = Disp[iDisp].GetNode();
		NodeDisplacement = Disp[iDisp].GetValue();

		for (int inode = 0; inode < nDispNode; inode++)
		{
			NodeIndex = Disp[iDisp].at(inode);
			nDispDof++;
			DispDof.at(NodeIndex * 2 + iDir) = nDispDof;
			IniDisplacement.at(NodeIndex * 2 + iDir) += NodeDisplacement.at(inode);
		}
		for (int igroup = 0; igroup < nGroup; igroup++)
		{
			int Type = Groups[igroup].GetType();
			int nEle = Groups[igroup].GetnElements();
			if (Type == 4)
			{
				Quadr *Elems=NULL;
				Elems = Groups[igroup].GetElement(*Elems);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					Elems[ielem].FillDof(DispDof);
				}
			}
		}
		nDispDof -= TotalDOF;
		LStiff.SetSize(TotalDOF, nDispDof);
		for (int igroup = 0; igroup < nGroup; igroup++)
		{
			int Type = Groups[igroup].GetType();
			int nEle = Groups[igroup].GetnElements();
			if (Type == 4)
			{
				Quadr *Elems=NULL;
				Elems = Groups[igroup].GetElement(*Elems);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					EStiff = Elems[ielem].GetStiff();
					EDof = Elems[ielem].GetDof();
					for (int idof = 0; idof < Type*nDof; idof++)
					{
						iedof = EDof.at(idof);
						for (int jdof = 0; jdof < Type*nDof; jdof++)
						{
							jedof = EDof.at(jdof);
							if (iedof != 0 && jedof != 0)
							{
								if (jedof>TotalDOF && iedof <= TotalDOF)
								{
									LStiff.at(iedof - 1, jedof - 1 - TotalDOF) += EStiff.at(idof, jdof);
								}
							}
						}
					}
				}
			}
		}
		for (int igroup = 0; igroup < nGroup; igroup++)
		{
			int Type = Groups[igroup].GetType();
			int nEle = Groups[igroup].GetnElements();
			if (Type == 4)
			{
				Quadr *Elems=NULL;
				Elems = Groups[igroup].GetElement(*Elems);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					Elems[ielem].FillDof(DegreeOfFreedom);
				}
			}
		}
		InitialDispLoad = InitialDispLoad - LStiff.Mult(NodeDisplacement);
	}
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems=NULL;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem].SetInitialDisplacement(IniDisplacement);
			}
		}
	}
}
void FemMain::GloableSolve()
{
	switch (ProblemType)
	{
	case Static:
		StaticSolve();
		break;
	case Model:
		ModelSolve();
		break;
	case DynamicStatic:
		DynamicStaticSolve();
		break;
	}
}
void FemMain::DynamicStaticSolve()
{
	Newmark Newmarker;
	Sor Sorer;
	CentralDifference CDiff;
	int iiter = 0;
	switch (SolveType)
	{
	case Explicit:
		CDiff.Init(dT, TotalDOF);
		for (int istep = 0; istep < nStep; istep++)
		{
			ComputeDOF();
			InitSolve();
			ComputeElementStiff();
			AssembleStiff();
			if (nInter > 0)
			{
				AssembleIStiff();
			}
			ApplyLoad();
			CountElement();
			switch (SolveMethod)
			{
			case MSor:
				Sorer.Init(Mass);
				break;
			case MLU:
				LUSolver.Decomposition(Mass);
				break;
			}
			TotalLoad = ExternalForce+ InitialStain + InteractLoad + InitialDispLoad;
			switch (SolveMethod)
			{
			case MSor:
				Sorer.Solve(TotalLoad, ResultSecond);
				break;
			case MLU:
				LUSolver.Solver(TotalLoad, ResultSecond);
				break;
			}
			//ResultSecond.Print();
			LResultZero= CDiff.Init(ResultSecond);

			EffictiveMass = CDiff.EffictiveMass(Mass, Damp);
			switch (SolveMethod)
			{
			case MSor:
				Sorer.Init(EffictiveMass);
				break;
			case MLU:
				LUSolver.Decomposition(EffictiveMass);
				break;
			}
			//EffictiveMass.Print();
			
			do
			{
				iiter++;

				//TotalLoad = ExternalForce + InitialStain + InteractLoad + InitialDispLoad;
				//EffictiveLoad = CDiff.EffictiveLoad(TotalLoad, Stiff, Mass, Damp,
				//	ResultZero, LResultZero);
				//switch (SolveMethod)
				//{
				//case MSor:
				//	Sorer.Solve(EffictiveLoad, ResultZero);
				//	break;
				//case MLU:
				//	LUSolver.Solver(EffictiveLoad, ResultZero);
				//	break;
				//}
				//CDiff.SolvePorcess(ResultZero, LResultZero, ResultFirst, ResultSecond);
				//ComputeElementStress();
				//CountElement();
				//SendResultToNode();

				if (nProces > 1)
				{
					
					InteractValue = ExchangeData();
					//for (int i = 0; i < nInter; i++)
					//{
					//	cout << "MyID=   " <<setw(10)<< Id <<"  Inter="<<setw(10)<<i<<"  "<< endl;
					//	InteractValue[i].Print();
					//}
					SetInteractResult(InteractValue);
				}
				TotalLoad = ExternalForce + InitialStain + InteractLoad + InitialDispLoad;
				EffictiveLoad = CDiff.EffictiveLoad(TotalLoad, Stiff, Mass, Damp,
					ResultZero, LResultZero);
				//cout << "ExternalForce    "<<setw(10)<<Id<<setw(10)<<iiter;
				//ExternalForce.Print();
				//cout << "InitialStain    " << setw(10) << Id << setw(10) << iiter;
				//InitialStain.Print();
				//cout << "InteractLoad    " << setw(10) << Id << setw(10) << iiter;
				//InteractLoad.Print();
				//cout << "InitialDispLoad    " << setw(10) << Id << setw(10) << iiter;
				//InitialDispLoad.Print();

				LResultZero = ResultZero;
				switch (SolveMethod)
				{
				case MSor:
					Sorer.Solve(EffictiveLoad, ResultZero);
					break;
				case MLU:
					LUSolver.Solver(EffictiveLoad, ResultZero);
					break;
				}
				CDiff.SolvePorcess(ResultZero, LResultZero, ResultFirst, ResultSecond);
				if (nProces > 1)
				{
					ConvergeCheck();
				}
				//cout << "ResultZero    " << setw(10) << Id << setw(10) << iiter;
				//ResultZero.Print();
				
				LResultFirst = ResultFirst;
				LResultSecond = ResultSecond;
				ComputeElementStress();
				
				SendResultToNode();
				GIDOutResult(iiter);
				
				//if (Id == 0)
				//{
					cout << setw(7) << "Iter=" << setw(7) << iiter
						<< setw(14) << "InterError=" <<setw(7)<< Error<<"  ";
					ShowTime();
				//}
			} while (Converge[1] == false && iiter < MaxIter);
		}
		break;
	case Implicit:
		Newmarker.IntSolver(dT);
		for (int istep = 0; istep < nStep; istep++)
		{
			ComputeDOF();
			InitSolve();
			ComputeElementStiff();
			AssembleStiff();
			ApplyLoad();
			TotalLoad = ExternalForce + InitialStain + InteractLoad + InitialDispLoad;

			EffictiveStiff = Newmarker.EffictiveStiff(Stiff, Mass, Damp);
			//EffictiveStiff.Print();
			switch (SolveMethod)
			{
			case MSor:
				Sorer.Init(EffictiveStiff);
				break;
			case MLU:
				LUSolver.Decomposition(EffictiveStiff);
				break;
			}
			do
			{
				iiter++;

				if (nProces > 1)
				{
					cout << "MyID=  " << Id;
				
					//InteractValue.Print();
					InteractValue = ExchangeData();
					//InteractValue.Print();
					SetInteractResult(InteractValue);
				}
				TotalLoad = ExternalForce + InitialStain + InteractLoad + InitialDispLoad;
				ShowTime();
				EffictiveLoad = Newmarker.EffictiveLoad(TotalLoad, LResultZero, LResultFirst, LResultSecond, Mass, Damp);
				//cout << "EffictiveLoad=";
				//EffictiveLoad.Print();

				switch (SolveMethod)
				{
				case MSor:
					Sorer.Solve(EffictiveLoad, ResultZero);
					break;
				case MLU:
					LUSolver.Solver(EffictiveLoad, ResultZero);
					break;
				}

				Newmarker.SolvePorcess(ResultZero, LResultZero, ResultFirst, LResultFirst, ResultSecond, LResultSecond);

				ConvergeCheck();
				cout << "Iter=    " << iiter << endl;
				LResultZero = ResultZero;
				LResultFirst = ResultFirst;
				LResultSecond = ResultSecond;
				ComputeElementStress();
				CountElement();
				SendResultToNode();
				GIDOutResult(iiter);
			} while (Converge[1] == false && iiter < MaxIter);
		}
		break;
	}
	
}
void FemMain::ModelSolve()
{
	ResultFirst.SetSize(TotalDOF);
	ResultZero.Set(1);
	LUSolve LUSolver;
	FloatArray Error(TotalDOF);
	double ErrorSum,ValueSum;
	LUSolver.Decomposition(Stiff);
	LUSolver.Inverse();
	for (int iStep = 0; iStep < nStep; iStep++)
	{
		LUSolver.Mult(ResultZero);
		Eigenvalues.at(iStep) = ResultZero.at(iStep);
		ResultZero = ResultZero.Times(1 / Eigenvalues.at(iStep));
		do
		{
			ResultFirst = ResultZero;
			LUSolver.Mult(ResultZero);
			Error = ResultZero - ResultFirst;
			ErrorSum = 0;
			ValueSum = 0;
			for (int i = 0; i < TotalDOF; i++)
			{
				ErrorSum += pow(Error.at(i), 2);
				ValueSum += pow(ResultZero.at(i), 2);
			}
			ErrorSum = ErrorSum / TotalDOF;
		} while (ErrorSum > Tolerance);
		//ResultZero.Print();
		ComputeElementStress();
		CountElement();
		SendResultToNode();
		GIDOutResult(iStep);
	}
	

}
void FemMain::StaticSolve()
{
	int iiter = 0;
	for (int istep = 0; istep < nStep; istep++)
	{
		ComputeDOF();
		InitSolve();
		ComputeElementStiff();
		AssembleStiff();
		ApplyLoad();

		TotalLoad = ExternalForce + InitialStain + InteractLoad + InitialDispLoad;
		ShowTime();
		switch (SolveMethod)
		{
		case MSor:
			LUSolver.Decomposition(Stiff);
			LUSolver.Solver(TotalLoad, ResultZero);
			LUSolver.Check(TotalLoad, ResultZero);
			break;
		case MLU:
			Soror.Init(Stiff);
			Soror.Solve(TotalLoad, ResultZero);
			break;
		}
		ComputeElementStress();
		CountElement();
		SendResultToNode();
		GIDOutResult(iiter);
		if (nProces>1)
		{
			do
			{
				iiter++;
				cout << "Myid= " << Id;
				//InteractValue.Print();
				InteractValue = ExchangeData();
				//InteractValue.Print();
				SetInteractResult(InteractValue);

				TotalLoad = ExternalForce + InitialStain + InteractLoad + InitialDispLoad;

				ShowTime();
				switch (SolveMethod)
				{
				case MSor:
					LUSolver.Decomposition(Stiff);
					LUSolver.Solver(TotalLoad, ResultZero);
					LUSolver.Check(TotalLoad, ResultZero);
					break;
				case MLU:
					Soror.Init(Stiff);
					Soror.Solve(TotalLoad, ResultZero);
					break;
				}
				ConvergeCheck();
				cout << "Iter=    " << iiter << endl;
				ComputeElementStress();
				CountElement();
				SendResultToNode();
				GIDOutResult(iiter);

			} while (Converge[1] == false && iiter < MaxIter);
		}
	}
	
}

bool *FemMain::ConvergeCheck()
{
	Error = 0;
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		//if (Id == 0)
		//{
			//InteractValueOld[iinter].Print();
		//cout << setw(10) << Id<<setw(10)<<iinter;
		//	InteractValue[iinter].Print();
		//}
		InteractValueOld[iinter] = InteractValue[iinter] - InteractValueOld[iinter];
		double Mean = InteractValue[iinter].Mean();
		if (Mean == 0)
		{
			Error += 1;
		}
		else if (Mean !=0)
		{
			Error += InteractValueOld[iinter].Norm() / Mean;
		}
		InteractValueOld[iinter] = InteractValue[iinter];
	}
	Converge[0] = false;
	Converge[1] = false;
	if (abs(Error) < Tolerance)
	{
		Converge[0] = true;
	}
	Converge[1] = Converge[0];
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		int AdjDomain = Inters[iinter].GetAdj();
		bool AdjConverge = false;
		MPI::COMM_WORLD.Send(&Converge[0], 1, MPI_C_BOOL, AdjDomain, TagCheck);
		MPI::COMM_WORLD.Recv(&AdjConverge, 1, MPI_C_BOOL, AdjDomain, TagCheck);
		Converge[1] = Converge[1] && AdjConverge;
	}

	//cout << "Error=   " << setw(20) << Error << setw(20) << "Converge[0]=   " << Converge[0] 
	//	<<"Converge[1]=   " << Converge[1] << endl;

	return Converge;
}

void FemMain::GIDOutResult(int istep)
{
	FloatArray Displacement, Stress, Strain;
	GiD_fBeginResult(FRes,"Displacement", "Static", istep, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL);
	for (int inode = 0; inode < nNode; inode++)
	{
		Displacement = Nodes[inode].GetDisplacement();
		GiD_fWrite2DVector(FRes, inode + 1, Displacement.at(0), Displacement.at(1));
	}
	GiD_fEndResult(FRes);
	GiD_fBeginResult(FRes,"Strain", "Static", istep, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL);
	for (int inode = 0; inode < nNode; inode++)
	{
		Strain = Nodes[inode].GetStrain();
		GiD_fWrite2DVector(FRes,inode + 1, Strain.at(0), Strain.at(1));
	}
	GiD_fEndResult(FRes);
	GiD_fBeginResult(FRes,"Stress", "Static", istep, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL);
	for (int inode = 0; inode < nNode; inode++)
	{
		Stress = Nodes[inode].GetStress();
		GiD_fWrite2DVector(FRes,inode + 1, Stress.at(0), Stress.at(1));
	}
	GiD_fEndResult(FRes);

}

void FemMain::CloseGidFile()
{
	GiD_fClosePostMeshFile(FMesh);
	GiD_fClosePostResultFile(FRes);
}

void FemMain::ComputeElementStress()
{

	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems=NULL;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem].SetResult(ResultZero);
				Elems[ielem].ComputeStress();
			}
		}
	}
}

void FemMain::CountElement()
{
	for (int inode = 0; inode < nNode; inode++)
	{
		Nodes[inode].ResetCount();
	}
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems = NULL;
			IntArray ENode(4);
			int NodeIndex;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ENode = Elems[ielem].GetNodeArray();
				for (int inode = 0; inode < Type; inode++)
				{
					NodeIndex = ENode.at(inode);
					Nodes[NodeIndex].AddCount();
				}
			}
		}
	}
}

void FemMain::SendResultToNode()
{
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems = NULL;
			IntArray ENode(4);
			FloatArray NodeStress(3);
			FloatArray NodeStrain(3);
			FloatArray Displacement(2);
			int NodeIndex;
			Elems = Groups[igroup].GetElement(*Elems);
			
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ENode = Elems[ielem].GetNodeArray();
				for (int inode = 0; inode < Type; inode++)
				{
					NodeIndex = ENode.at(inode);
					NodeStrain = Elems[ielem].GetStrain(inode);
					NodeStress = Elems[ielem].GetStress(inode);
					Displacement = Elems[ielem].GetDisplacement(inode);
					//cout << "NodeStrain" << setw(10) << NodeIndex;
					//NodeStrain.Print();
					//cout << "NodeStress"<<setw(10)<<NodeIndex; 
					//NodeStress.Print();
					Nodes[NodeIndex].SetStrain(NodeStrain);
					Nodes[NodeIndex].SetStress(NodeStress); 
					Nodes[NodeIndex].SetDisplacement(Displacement);
				}
			}
		}
	}
}

IntArray FemMain::GetInteractNode(int iinter)
{
	return Inters[iinter].GetRemote();
}
FloatArray FemMain::GetInteractResult(IntArray & InteractNode)
{
	FloatArray InteractResult;
	FloatArray NodeDisplacement;
	int nInterNode = InteractNode.GetSize();
	InteractResult.SetSize(nInterNode*nDof);
	for (int inode = 0; inode < nInterNode; inode++)
	{
		int NodeIdx = InteractNode.at(inode);
		NodeDisplacement = Nodes[NodeIdx].GetDisplacement();
		for (int iDof = 0; iDof < nDof; iDof++)
		{
			InteractResult.at(inode * 2 + iDof) = NodeDisplacement.at(iDof);
		}
	}
	return InteractResult;
}
void FemMain::AssembleIStiff()
{
	IntArray LocalNode, RemoteNode;

	FloatMatrix EStiff;
	int iedof, jedof;
	IntArray EDof;
	int NodeIndex, nInterDof, nInterNode;

	InterDof = DegreeOfFreedom;
	nInterDof = TotalDOF;
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		LocalNode = Inters[iinter].GetLocal();
		RemoteNode = Inters[iinter].GetRemote();
		nInterNode = LocalNode.GetSize();
		for (int inode = 0; inode < nInterNode; inode++)
		{
			NodeIndex = LocalNode.at(inode);
			for (int iDof = 0; iDof < nDof; iDof++)
			{
				nInterDof++;
				InterDof.at(NodeIndex * 2 + iDof) = nInterDof;
			}
		}
		//cout << setw(20) << "InteractResult[iinter]=" << setw(10) << "ID" << setw(10) << Id
		//	<<setw(10)<<iinter;
		//InteractResult[iinter].Print();
	}

	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems = NULL;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem].FillDof(InterDof);
			}
		}
	}
	nInterDof -= TotalDOF;

	IStiff.SetSize(TotalDOF, nInterDof);
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems = NULL;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				EStiff = Elems[ielem].GetStiff();
				EDof = Elems[ielem].GetDof();
				for (int idof = 0; idof < Type*nDof; idof++)
				{
					iedof = EDof.at(idof);
					for (int jdof = 0; jdof < Type*nDof; jdof++)
					{
						jedof = EDof.at(jdof);
						if (iedof != 0 && jedof != 0)
						{
							if (jedof>TotalDOF && iedof <= TotalDOF)
							{
								IStiff.at(iedof - 1, jedof - 1 - TotalDOF) += EStiff.at(idof, jdof);
							}
						}
					}
				}
			}
		}
	}
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems = NULL;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem].FillDof(DegreeOfFreedom);
			}
		}
	}
}
void FemMain::SetInteractResult(FloatArray * InteractResult)
{
	IntArray LocalNode, RemoteNode;

	int NodeIndex, nInterDof, nInterNode;

	InterDisplace.Clear();
	InteractLoad.Clear();
	
	nInterDof = TotalDOF;
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		LocalNode = Inters[iinter].GetLocal();
		RemoteNode = Inters[iinter].GetRemote();
		nInterNode = LocalNode.GetSize();
		for (int inode = 0; inode < nInterNode; inode++)
		{
			NodeIndex = LocalNode.at(inode);
			for (int iDof = 0; iDof < nDof; iDof++)
			{
				InterDisplace.at(NodeIndex * 2 + iDof) += InteractResult[iinter].at(inode * 2 + iDof);
			}
		}
		//cout << setw(20) << "InteractResult[iinter]=" << setw(10) << "ID" << setw(10) << Id
		//	<<setw(10)<<iinter;
		//InteractResult[iinter].Print();
	}
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr *Elems = NULL;
			Elems = Groups[igroup].GetElement(*Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem].SetInitialDisplacement(InterDisplace);
			}
		}
	}
	FloatArray iInteractResult(0);
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		iInteractResult =iInteractResult.Cat(InteractResult[iinter]);
		//if (Id == 1)
		//{
		//	cout << "InteractResult[iinter]" << setw(10) << Id;
		//	InteractResult[iinter].Print();
		//}
	}
	//if (Id == 1)
	//{
	//	cout << "iInteractResult" << setw(10) << Id;
	//	iInteractResult.Print();
	//}
	InteractLoad = IStiff.Mult(iInteractResult);
	//if (Id == 1)
	//{
	//	//cout << "IStiff" << setw(10) << Id;
	//	//IStiff.Print();
	//	cout << "InteractLoad" << setw(10) << Id;
	//	InteractLoad.Print();
	//}
}

void FemMain::GetSize(int &NProces)
{
	nProces = NProces;
}

void FemMain::GetID(int &MyID)
{
	Id = MyID;
}

FloatArray * FemMain::ExchangeData()
{
	int *RemoteNode;
	double *Value;
	MPI::COMM_WORLD.Barrier();
	for (int iinter = 0; iinter < nInter; iinter++)
	{
		InteractNode = GetInteractNode(iinter);
		int AdjDomain = Inters[iinter].GetAdj();
		int size = InteractNode.GetSize();

		RemoteNode = new int[size]();
		Value = new double[size*nDof]();
		RemoteNode = InteractNode.GetValue();
		//cout << AdjDomain << endl;
		MPI::COMM_WORLD.Send(RemoteNode, size, MPI_INTEGER, AdjDomain, TagNode);
		MPI::COMM_WORLD.Recv(RemoteNode, size, MPI_INTEGER, AdjDomain, TagNode);
		for (int i = 0; i < size; i++)
		{
			InteractNode.at(i) = RemoteNode[i];
		}
		//InteractNode.Print();
		InteractValue[iinter] = GetInteractResult(InteractNode);
		//cout << setw(10) << "ID="<<setw(10)<<Id<<setw(10)<<"before"<<iinter;
		//InteractValue[iinter].Print();
		Value = InteractValue[iinter].GetValue();
		MPI::COMM_WORLD.Send(Value, size*nDof, MPI_DOUBLE, AdjDomain, TagValue);
		MPI::COMM_WORLD.Recv(Value, size*nDof, MPI_DOUBLE, AdjDomain, TagValue);
		for (int i = 0; i < size*nDof; i++)
		{
			InteractValue[iinter].at(i) = Value[i];
		}
		//cout << setw(10) << "ID=" << setw(10) << Id << setw(10) << iinter;
		//InteractValue[iinter].Print();
		delete[]RemoteNode, Value;
	}
	return InteractValue;
}