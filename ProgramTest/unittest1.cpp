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
			A.Print();
			IntArray *B;
			B = new IntArray(10);
			B->Print();
			*B = A.Append(*B);
			*B = B->Add(A);
			*B = A.Append(*B);
			B->Print();
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

	};
}