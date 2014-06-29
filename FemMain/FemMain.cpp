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
	FemMain Fem;
	string text1,text2;
	stringstream stream;
	ifstream inp("1.inp");
	IntArray InteractNode;
	FloatArray InteractValue,InteractValueOld;
	bool *Converge;
	double Error;
	int MaxIter;
	int NProces, MyID;
	MPI::Init(argc, argv);
	NProces = MPI::COMM_WORLD.Get_size();
	MyID = MPI::COMM_WORLD.Get_rank();
	Fem.GetSize(NProces);
	Fem.GetID(MyID);
	getline(inp, text1);
	getline(inp, text1);
	getline(inp, text2);
	if (MyID == 0)
	{
		Fem.WorkDir() = text1;
	}
	else if (MyID == 1)
	{
		Fem.WorkDir() = text2;
	}
	
	inp.close();
	MPI::COMM_WORLD.Barrier();
	GiD_PostInit();

	Fem.ShowTime();

	Fem.ReadFiles();

	Fem.OpenGidFile(); 


	Fem.GIDOutMesh();

	int nStep = Fem.GetStep();
	int iiter = 0;
	MaxIter = Fem.GetMaxIter();
	Fem.ShowTime();
	MPI::COMM_WORLD.Barrier();

	for (int istep = 0; istep < nStep; istep++)
	{
		Fem.ComputeDOF();

		Fem.InitSolve();

		Fem.ComputeElementStiff();


		Fem.AssembleStiff();
	

		Fem.ApplyLoad();

		
		Fem.Solve();

		Fem.ComputeElementStress();


		Fem.CountElement();


		Fem.SendResultToNode();

		Fem.GIDOutResult(iiter);
		Converge = new bool[2]();
		Converge[0] = Converge[1] = false;
		do 
		{
			iiter++;
			
			InteractValue=Fem.ExchangeData();
			Fem.SetInteractResult(InteractValue);
			if (Converge[0] == false)
			{
				Fem.Solve();
			}

			Converge = Fem.ConvergeCheck();

			cout << "Iter=    " << iiter << endl;

			Fem.ComputeElementStress();

			Fem.CountElement();

			Fem.SendResultToNode();

			Fem.GIDOutResult(iiter);


		} while (Converge[1] == false && iiter<MaxIter);

	}
	Fem.CloseGidFile();

	GiD_PostDone();
	MPI::COMM_WORLD.Barrier();
	MPI::Finalize();
	return 0;
}
void FemMain::ShowTime()
{
	time(&tNow);
	localtime_s(&tmLocal, &tNow);
	cout << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon + 1 << "-" <<
		tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0') << setw(2) << tmLocal.tm_min << ":"
		<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
	//chk << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon << "-" <<
	//	tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0') << setw(2) << tmLocal.tm_min << ":"
	//	<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
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
	stream >> nDim >> nNode >> nGroup >> nElem >> nMat >> nStep >> nDof >> MaxIter >> Tolerance;

	cout << Tolerance << endl;

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
		Mat--;
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
		Nodes[inode].Init(Idx, *Coors);
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
		int Dof = Groups[igroup].GetDof();
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
				Quadrs[ielem].Init(4, Mat, Idx, igroup,Dof, &Nodes);
			}
			for (int ielem = 0; ielem < nElem; ielem++)
			{
				Quadrs[ielem].Print();
			}
			Groups[igroup].FillElement(Quadrs);
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
	stream >> nTotalLoad >> nFace >> nVolumn >> nConcentrate;
	Faces = new Face[nFace]();
	Cons = new Concentrate[nConcentrate]();
	Vols = new Volumn[nVolumn]();
	for (int iface = 0; iface < nFace; iface++)
	{
		int nEle = 0, nNode = 0, Dir = 0,AdjElem=0;
		double StartC, EndC, StartV, EndV;
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);
		loa.getline(str, MAXCHAR);

		stream.str("");
		stream.clear();
		stream << str;
		stream >> nEle >> nNode >> Dir >> StartC >> EndC >> StartV >> EndV;
		Dir--;
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
				Lines[ielem].Init(2, NULL, ielem, NULL,NULL, &Nodes);
				Lines[ielem].AtAdjElem() = AdjElem;
			}
			Faces[iface].Set(Lines);
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
		Dir--;
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
		Vols[ivol].Init(ivol, nGroup, *group, Acc, Dir);
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
		Dir--;
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
		Cons[icon].Init(icon, nNode,* Node, *Value, Dir);
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
		Dir--;
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
		Fix[ifix].Init(ifix, nNode,* Nodes, Dir);
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
		Dir--;
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
		Disp[idisp].Init(nNode, Dir, *Node, *Value);
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
		AdjDomain--;
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
		Inters[iinter].Init(nNode, AdjDomain,* Local,* Remote);
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
	FloatArray *Coor;
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int type,nEle;
		int *nid;
		type = Groups[igroup].GetType();
		nEle = Groups[igroup].GetnElements();
		GiD_fBeginMesh(FMesh, GroupName[type - 1].c_str(), GiD_Dimension(nDim), GiD_ElementType(type), type);
		if (igroup == 0)
		{
			GiD_fBeginCoordinates(FMesh);
			Coor = new FloatArray(nDim);
			for (int inode = 0; inode < nNode; inode++)
			{
				Coor = &Nodes[inode].GetCoordinate();
				if (nDim == 2)
				{
					GiD_fWriteCoordinates2D(FMesh,inode + 1, Coor->at(0), Coor->at(1));
				}
				else
				{
					GiD_fWriteCoordinates(FMesh,inode + 1, Coor->at(0), Coor->at(1), Coor->at(2));
				}
			}
			GiD_fEndCoordinates(FMesh);
		}
		GiD_fBeginElements(FMesh);
		if (type == 4)
		{
			Quadr **Elems;
			IntArray Nodes;
			Elems = new Quadr*[nEle];
			nid = new int[type];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int iele = 0; iele < nEle; iele++)
			{
				Nodes = Elems[iele]->GetNodeArray();
				for (int inode = 0; inode < type; inode++)
				{
					nid[inode] = Nodes.at(inode)+1;
				}
				GiD_fWriteElement(FMesh,Elems[iele]->GetIndex(), nid);
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
	cout << "DegreeOfFreedom      " ;
	DegreeOfFreedom.Print();
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr **Elems;
			Elems = new Quadr*[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem]->FillDof(DegreeOfFreedom);
			}
		}
	}
}

void FemMain::InitSolve()
{
	Stiff.SetSize(TotalDOF, TotalDOF);
	ResultZero .SetSize(TotalDOF);
	ExternalForce .SetSize(TotalDOF);
	InitialStain .SetSize(TotalDOF);
	InitialDispLoad .SetSize(TotalDOF);
	IniDisplacement.SetSize(nNode*nDof);
	InterDisplace.SetSize(nNode*nDof);
	TotalLoad.SetSize(TotalDOF);
	InteractLoad.SetSize(TotalDOF);
	Converge = new bool[2];
	int nRow, nCol, NonZero;
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
			Quadr **Elems;
			Elems = new Quadr*[nEle];
			EDof = new IntArray(Type*Dof);
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				*EDof=Elems[ielem]->GetDof();
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
			Quadr **Elems;
			IntArray ENode(Type);
			FloatMatrix *Coor;
			Material *Mat;
			Elems = new Quadr*[nEle];
			Coor = new FloatMatrix (Type,nDim);
			Mat = new Material(Mats[Groups[igroup].GetMaterial()]);
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ENode = Elems[ielem]->GetNodeArray();
				Elems[ielem]->SetMaterial(*Mat);
				for (int inode = 0; inode < Type; inode++)
				{
					for (int idim = 0; idim < nDim; idim++)
					{
						Coor->at(inode, idim) = Nodes[ENode.at(inode)].GetCoordinate().at(idim);
					}
				}
				Elems[ielem]->SetCoor(*Coor);
				Elems[ielem]->ComputeStiff();
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
			Quadr **Elems;
			FloatMatrix EStiff(Type*nDof, Type*nDof);
			Elems = new Quadr*[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			IntArray Dof(Type*nDof);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Dof = Elems[ielem]->GetDof();
				EStiff = Elems[ielem]->GetStiff();
				for (int idof = 0; idof < Type*nDof; idof++)
				{
					iedof = Dof.at(idof);
					for (int jdof = 0; jdof < Type*nDof; jdof++)
					{
						jedof = Dof.at(jdof);
						if (iedof != 0 && jedof != 0)
						{
							Stiff.at(iedof - 1, jedof - 1) += EStiff.at(idof, jdof);
						}
					}
				}
			}

		}
	}
	//Stiff.Print();
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
			Line *Lines;
			int Type = 2;
			FloatMatrix *Coor;
			Lines = new Line[nEle];
			IntArray *ENode;
			ENode = new IntArray(2);
			double Length;
			FloatArray *PVal, *Normal, *Shape, *GaussCoor, *TELoad, *ELoad;
			PVal = new FloatArray(2);
			Normal = new FloatArray(3);
			Shape = new FloatArray(2);
			ELoad = new FloatArray(2);
			TELoad = new FloatArray(2);
			GaussCoor = new FloatArray(1);
			GaussPoint *B;
			B = new GaussPoint();
			double Det,PLoad; 
			Lines = Faces[iface].GetLines();
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ELoad->Clear();
				*ENode = Lines[ielem].GetNodeArray();
				Coor = new FloatMatrix(Type, nDim);
				for (int inode = 0; inode < Type; inode++)
				{
					for (int idim = 0; idim < nDim; idim++)
					{
						Coor->at(inode, idim) = Nodes[ENode->at(inode)].GetCoordinate().at(idim);
					}
					PVal->at(inode) = (Coor->at(inode, Dir) - StartC) / (EndC - StartC)*(EndV - StartV) + StartV;
				}
				Length = sqrt(pow(Coor->at(0, 0) - Coor->at(1, 0), 2) +
					pow(Coor->at(0, 1) - Coor->at(1, 1), 2));
				Normal->at(0) = (Coor->at(0, 0) - Coor->at(1, 0)) / Length;
				Normal->at(1) = (Coor->at(0, 1) - Coor->at(1, 1)) / Length;
				FloatArray Descartes(3);
				Descartes.at(2) = 1;
				*Normal = Normal->Cross(Descartes);
				double ksi, weight;
				for (int iksi = 0; iksi < 2; iksi++)
				{
					ksi = Gauss2[iksi];
					weight = Weight2[iksi];
					Shape->at(0) = 0.5*(1 - ksi);
					Shape->at(1) = 0.5*(1 + ksi);
					Det = 0.5*Length;
					PLoad = PVal->Dot(Shape);
					PLoad = PLoad*Det*weight;
					TELoad = &Shape->Times(PLoad);
					ELoad->at(0) += TELoad->at(0);
					ELoad->at(1) += TELoad->at(1);
				}
				for (int inode = 0; inode < Type; inode++)
				{
					int NodeIdx = ENode->at(inode);
					for (int idof = 0; idof < nDof; idof++)
					{
						if (DegreeOfFreedom.at(NodeIdx * 2 + idof) != 0)
						{
							ExternalForce.at(DegreeOfFreedom.at(NodeIdx * 2 + idof) - 1) -= ELoad->at(idof)*Normal->at(idof);
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
		Quadr **Elems;
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
				Elems = new Quadr*[nEle];
				Elems = Groups[GroupIdx].GetElement(**Elems);
				FloatArray *Eload = new FloatArray(2);
				FloatArray GaussT(2);
				FloatArray Shape(Type);
				FloatArray TLoad(Type);
				IntArray Nodes(Type);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					TLoad.Clear();
					Nodes = Elems[ielem]->GetNodeArray();
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
							Elems[ielem]->ComputeJacobi(B);
							Shape = Elems[ielem]->GetShape();
							Det = Elems[ielem]->GetDet();
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
				Quadr **Elems;
				Elems = new Quadr*[nEle];
				Elems = Groups[igroup].GetElement(**Elems);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					Elems[ielem]->FillDof(DispDof);
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
				Quadr **Elems;
				Elems = new Quadr*[nEle];
				Elems = Groups[igroup].GetElement(**Elems);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					EStiff = Elems[ielem]->GetStiff();
					EDof = Elems[ielem]->GetDof();
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
				Quadr **Elems;
				Elems = new Quadr*[nEle];
				Elems = Groups[igroup].GetElement(**Elems);
				for (int ielem = 0; ielem < nEle; ielem++)
				{
					Elems[ielem]->FillDof(DegreeOfFreedom);
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
			Quadr **Elems;
			Elems = new Quadr *[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem]->SetInitialDisplacement(IniDisplacement);
			}
		}
	}
}

void FemMain::Solve()
{
	TotalLoad = ExternalForce + InitialStain + InteractLoad + InitialDispLoad;
	//cout << "ExternForce";
	//ExternalForce.Print();
	//cout << "InteractLoad";
	//InteractLoad.Print();
	ShowTime();
	LUSolver = new LUSolve();
	LUSolver->Decomposition(Stiff);
	LUSolver->Solver(TotalLoad, ResultZero);
	LUSolver->Check(TotalLoad, ResultZero);
}

bool *FemMain::ConvergeCheck()
{
	InteractValueOld = InteractValue-InteractValueOld;
	Error = InteractValueOld.Norm()/InteractValue.Mean();
	Converge[0] = false;
	Converge[1] = false;
	if (abs(Error) < Tolerance)
	{
		Converge[0] = true;
	}
	int AdjDomain = Inters[0].GetAdj();
	
	bool AdjConverge = false;
	MPI::COMM_WORLD.Send(&Converge[0], 1, MPI_C_BOOL, AdjDomain, TagCheck);
	MPI::COMM_WORLD.Recv(&AdjConverge, 1, MPI_C_BOOL, AdjDomain, TagCheck);
	
	Converge[1] = Converge[0] && AdjConverge;

	InteractValueOld = InteractValue;

	cout << "Error=   " << setw(20) << Error << setw(20) << "Converge[0]=   " << Converge[0] 
		<<"Converge[1]=   " << Converge[1] << endl;

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
			Quadr **Elems;
			Elems = new Quadr *[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem]->SetResult(ResultZero);
				Elems[ielem]->ComputeStress();
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
			Quadr **Elems;
			IntArray ENode(4);
			int NodeIndex;
			Elems = new Quadr *[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ENode = Elems[ielem]->GetNodeArray();
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
			Quadr **Elems;
			IntArray ENode(4);
			FloatArray NodeStress(3);
			FloatArray NodeStrain(3);
			FloatArray Displacement(2);
			int NodeIndex;
			Elems = new Quadr *[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				ENode = Elems[ielem]->GetNodeArray();
				for (int inode = 0; inode < Type; inode++)
				{
					NodeIndex = ENode.at(inode);
					NodeStrain = Elems[ielem]->GetStrain(inode);
					NodeStress = Elems[ielem]->GetStress(inode);
					Displacement = Elems[ielem]->GetDisplacement(inode);
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

IntArray FemMain::GetInteractNode()
{
	return Inters[0].GetRemote();
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
void FemMain::SetInteractResult(FloatArray & InteractResult)
{
	IntArray LocalNode, RemoteNode;
	LocalNode = Inters[0].GetLocal();
	RemoteNode = Inters[0].GetRemote();

	FloatMatrix A;
	FloatMatrix IStiff, EStiff;
	int iedof, jedof;
	IntArray EDof;
	int NodeIndex, nInterDof, nInterNode;
	IntArray InterDof;

	InterDisplace.Clear();
	InteractLoad.Clear();

	nInterNode = LocalNode.GetSize();
	A.SetSize(TotalDOF, TotalDOF);
	InterDof = DegreeOfFreedom;
	nInterDof = TotalDOF;
	for (int inode = 0; inode < nInterNode; inode++)
	{
		NodeIndex = LocalNode.at(inode);
		for (int iDof = 0; iDof < nDof; iDof++)
		{
			nInterDof++;
			InterDof.at(NodeIndex * 2 + iDof) = nInterDof;
			InterDisplace.at(NodeIndex * 2 + iDof) += InteractResult.at(inode * 2 + iDof);
		}
	}
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr **Elems;
			Elems = new Quadr*[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem]->FillDof(InterDof);
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
			Quadr **Elems;
			Elems = new Quadr*[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				EStiff = Elems[ielem]->GetStiff();
				EDof = Elems[ielem]->GetDof();
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
			Quadr **Elems;
			Elems = new Quadr*[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem]->FillDof(DegreeOfFreedom);
			}
		}
	}
	//cout << "InterDisplace   ";
	//InterDisplace.Print(); 
	for (int igroup = 0; igroup < nGroup; igroup++)
	{
		int Type = Groups[igroup].GetType();
		int nEle = Groups[igroup].GetnElements();
		if (Type == 4)
		{
			Quadr **Elems;
			Elems = new Quadr *[nEle];
			Elems = Groups[igroup].GetElement(**Elems);
			for (int ielem = 0; ielem < nEle; ielem++)
			{
				Elems[ielem]->SetInitialDisplacement(InterDisplace);
			}
		}
	}
	InteractLoad = InteractLoad - IStiff.Mult(InteractResult);
}

void FemMain::GetSize(int &NProces)
{
	nProces = NProces;
}

void FemMain::GetID(int &MyID)
{
	Id = MyID;
}

FloatArray FemMain::ExchangeData()
{
	int *RemoteNode;
	double *Value;
	MPI::COMM_WORLD.Barrier();
	InteractNode = GetInteractNode();
	
	int AdjDomain = Inters[0].GetAdj();
	int size = InteractNode.GetSize();
	InteractValue.SetSize(size*nDof);

	RemoteNode = new int[size]();
	Value = new double[size*nDof]();
	RemoteNode = InteractNode.GetValue();
	cout << AdjDomain << endl;
	MPI::COMM_WORLD.Send(InteractNode.GetValue(), size, MPI_INTEGER, AdjDomain, TagNode);
	MPI::COMM_WORLD.Recv(RemoteNode, size, MPI_INTEGER, AdjDomain, TagNode);
	for (int i = 0; i < size; i++)
	{
		InteractNode.at(i) = RemoteNode[i];
	}

	InteractValue = GetInteractResult(InteractNode);

	Value = InteractValue.GetValue();
	MPI::COMM_WORLD.Send(InteractValue.GetValue(), size*nDof, MPI_DOUBLE, AdjDomain, TagValue);
	MPI::COMM_WORLD.Recv(Value, size*nDof, MPI_DOUBLE, AdjDomain, TagValue);
	for (int i = 0; i < size*nDof; i++)
	{
		InteractValue.at(i) = Value[i];
	}
	return InteractValue;
}