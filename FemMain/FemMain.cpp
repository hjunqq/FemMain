#include "Array.h"
int main()
{
	FloatArray A(10);
	for (int i = 0; i < 10; i++)
	{
		A.at(i) = 1;
	}
	A.Print();
	FloatMatrix B(10, 10);
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			B.at(i, j) = i + j;
		}
	}
	B.Print();
	B = B.Trans();
	B.Print();
	A = B.Mult(&A);
	A.Print();
	B = B.Mult(&B);
	B.Print();
	return 0;
}