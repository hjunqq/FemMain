#include "stdafx.h"
#include "CppUnitTest.h"
#include "../FemMain/Array.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace ProgramTest
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			// TODO:  在此输入测试代码
			IntArray A(10);
			for (int i = 0; i < 10; i++)
			{
				A.at(i) = i;
			}
			IntArray *B;
			B = new IntArray(10);
			*B = A.Append(*B);
			*B = B->Add(A);
			*B = A.Append(*B);
			Assert::AreEqual(30, B->GetSize());
		}
		TEST_METHOD(TestMethod2)
		{
			FloatArray A(10);
			for (int i = 0; i < 10; i++)
			{
				A.at(i) = i;
			}
			FloatMatrix B(10, 10);
			for (int i = 0; i < 10; i++)
			{
				for (int j = 0; j < 10; j++)
				{
					B.at(i, j) = i + j;
				}
			}
			B = B.Trans();
			A = B.Mult(&A);
		}
		TEST_METHOD(TestMethod3)
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
	};
}