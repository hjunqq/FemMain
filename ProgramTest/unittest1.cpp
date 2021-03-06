#include "stdafx.h"
#include "CppUnitTest.h"
#include "../FemMain/Array.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace ProgramTest
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestIntArray)
		{
			// TODO:  �ڴ�������Դ���
			IntArray A(10);
			for (int i = 0; i < 10; i++)
			{
				A.at(i) = i;
			}
			IntArray B;
			B.SetSize(10);
			B = A;
			Assert::AreEqual(B.at(5), A.at(5));
			B = B + A;
			Assert::AreEqual(10, B.at(5));
			B = B - A;
			Assert::AreEqual(5, B.at(5));
		}
		TEST_METHOD(TestFloatArray)
		{
			FloatArray A(10);
			for (int i = 0; i < 10; i++)
			{
				A.at(i) = i;
			}
			A = A.Times(2);
			Assert::AreEqual(2.0, A.at(1));
			double sum = A.Dot(&A);
			FloatMatrix B(10, 10);
			for (int i = 0; i < 10; i++)
			{
				for (int j = 0; j < 10; j++)
				{
					B.at(i, j) = i*10 + j;
				}
			}
			B = B.Trans();
			A = B.Mult(A);
		}
		TEST_METHOD(TestCSR)
		{
			
			double Values;
			int RowIdx[7] = { 0, 3, 5, 8, 12, 14, 15 };
			int ColIdx[15] = { 0, 3, 5, 0, 1, 1, 2, 5, 0, 3, 4, 5, 1, 4, 5 };
			int nRow, nCol,NonZero;
			nRow = 6;
			nCol = 6;
			NonZero = 15;
			CSRMatrix A;
			A.Init(NonZero, nRow, nCol, RowIdx, ColIdx);
			Values = 0;
			for (int iRow = 0; iRow < nRow; iRow++)
			{
				int i = iRow;
				for (int iCol = RowIdx[iRow]; iCol < RowIdx[iRow + 1]; iCol++)
				{
					int j = ColIdx[iCol];
					Values++;
					A.at(i, j) = Values;
				}				
			}
			Assert::AreEqual(13.0, A.at(4, 1));
		}
		TEST_METHOD(TestMatrix)
		{
			IntMatrix A(2,2);
			A.at(0, 0) = 1;
			A.at(0, 1) = 2;
			A.at(1, 0) = 3;
			A.at(1, 1) = 4;
			double det=A.Determinant();
			Assert::AreEqual(-2.0, det);
		}
		/*TEST_METHOD(TestNode)
		{
			Node A;
			FloatArray *Coor;
			Coor = new FloatArray(3);
			for (int icoor = 0; icoor < 3; icoor++)
			{
				Coor->at(icoor) = icoor;
			}
			A.Init(0, Coor);
			A.Print();
		}
		TEST_METHOD(TestGroup)
		{
			Group A;
			IntArray *Elems;
			Elems = new IntArray(10);
			for (int i = 0; i < 10; i++)
			{
				Elems->at(i) = i + 1;
			}
			A.Init(0, 10, 1,1,2);
			A.FillElement(Elems);

			Quadr C;
			IntArray B(4);
			C.Print();
			C.Init(4, 9, 1, 1, Elems);
			Assert::AreEqual(3, C.GetNode(2));
		}
		TEST_METHOD(TestElement)
		{
			Line *A;
			A = new Line();
			IntArray *Nodes;
			Nodes = new IntArray(2);
			for (int i = 0; i < 10; i++)
			{
				Nodes->at(i) = i;
			}
			A->Init(1, 1, 1, 1, Nodes);
			Line B;
			B = *A;
			Assert::AreEqual(1, B.GetNode(2));*/
		//}
	};
	TEST_CLASS(TimeTest)
	{
	public:
		TEST_METHOD(TimeAdd)
		{
			Time coding(2, 40);
			Time fixing(5, 55);
			Time total;

			total = fixing;
			total.Show();
		}
	};
}