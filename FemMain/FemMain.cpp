#include "Array.h"
int main()
{
	double Values;
	int RowIdx[7] = { 0, 3, 5, 8, 12, 14, 15 };
	int ColIdx[15] = { 0, 3, 5, 0, 1, 1, 2, 5, 0, 3, 4, 5, 1, 4, 5 };
	int nRow, nCol, NonZero;
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
	A.Print();
	return 0;
}