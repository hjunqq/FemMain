#include "Array.h"

IntArray::IntArray()
{
	size = 0;
}


IntArray::~IntArray()
{
}


// 构造大小为n的数组
IntArray::IntArray(int n)
{
	size = n;
	Values = new int[n]();
}


 //返回第i个函数值	
int & IntArray::at(int i)
{
	return Values[i];
}



// 两个向量直接叠加
IntArray IntArray::Append(IntArray & TArray)
{
	IntArray Temp;
	int nNew;
	nNew = size + TArray.GetSize();
	Temp.Values = new int[nNew]();
	Temp.size = nNew;
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = Values[i];
	}
	for (int i = size; i < nNew; i++)
	{
		Temp.Values[i] = TArray.at(i - size);
	}
	return Temp;
}


// 数组相加
IntArray IntArray::Add(IntArray & TArray)
{
	IntArray Temp;
	Temp.Values = new int[size]();
	Temp.size =size;
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i]=Values[i] + TArray.at(i);
	}
	return Temp;
}


// 获取数组大小
int IntArray::GetSize()
{
	return size;
}


// 打印变量内容
void IntArray::Print()
{
	cout << "The Values are:" << endl;
	for (int i = 0; i < size; i++)
	{
		cout << setw(10) << Values[i];
	}
	cout << endl;
}


FloatArray::FloatArray()
{
	size = 0;
}


FloatArray::~FloatArray()
{
}


// 返回第i个元素
double & FloatArray::at(int i)
{
	return Values[i];
}


// 构造大小为n的数组
FloatArray::FloatArray(int n)
{
	size = n;
	Values = new double[n]();
}


// 返回数组大小
int FloatArray::GetSize()
{
	return size;
}


// 打印数组
void FloatArray::Print()
{
	cout << "FloatArray" << endl;
	for (int i = 0; i < size; i++)
	{
		cout << setw(15) << Values[i] << endl;
	}
	cout << endl;
}


// 乘以标量
FloatArray FloatArray::Times(double Scalar)
{
	FloatArray TArray;
	TArray.size = size;
	for (int i = 0; i < size; i++)
	{
		TArray.at(i) = Values[i] * Scalar;
	}
	return TArray;
}


// 点乘
double FloatArray::Dot(FloatArray * B)
{
	double sum=0.0;
	for (int i = 0; i < size; i++)
	{
		sum += Values[i] * B->at(i);
	}
	return sum;
}


// 清空
int FloatArray::Clear()
{
	delete[] Values;
	size = 0;
	return 0;
}


// 数组相加
FloatArray FloatArray::Plus(FloatArray * B)
{
	FloatArray Temp;
	Temp.size = size;
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = Values[i] + B->at(i);
	}
	return Temp;
}


// 数组相减
FloatArray FloatArray::Minus(FloatArray * B)
{
	FloatArray Temp;
	Temp.size = size;
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = Values[i] - B->at(i);
	}
	return Temp;
}


FloatMatrix::FloatMatrix()
{
	m = 0;
	n = 0;
}


FloatMatrix::~FloatMatrix()
{
}


IntMatrix::IntMatrix()
{
	m = 0;
	n = 0;
}


IntMatrix::~IntMatrix()
{
}


// 构造矩阵
IntMatrix::IntMatrix(int i, int j)
{
	m = i; n = j;
	Values = new int[m*n]();
}


// 返回第m行n列元素
int & IntMatrix::at(int i, int j)
{
	return Values[i*n + j];
}


// 构造m行n列的矩阵
FloatMatrix::FloatMatrix(int i, int j)
{
	m = i; n = j;
	Values = new double[m*n]();
}


// 返回第m行n列的值
double & FloatMatrix::at(int i, int j)
{
	return Values[i*n + j];
}


// 转置
FloatMatrix FloatMatrix::Trans()
{
	FloatMatrix Temp(n,m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i, j) = this->at(j, i);
		}
	}
	return Temp;
}


// 乘以数组
FloatArray FloatMatrix::Mult(FloatArray * B)
{
	FloatArray Temp(m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i) += B->at(j)*this->at(i,j);
		}
	}
	return Temp;
}

 // 乘以矩阵
FloatMatrix FloatMatrix::Mult(FloatMatrix * B)
{
	FloatMatrix Temp(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < B->n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				Temp.at(i, j) = this->at(i, k)*B->at(k,j);
			}			
		}
	}
	return Temp;
}

 // 矩阵相加
FloatMatrix FloatMatrix::Plus(FloatMatrix * B)
{
	FloatMatrix Temp(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i, j) = this->at(i, j) + B->at(i, j);
		}
	}
	return Temp;
}

 // 矩阵相减
FloatMatrix FloatMatrix::Minus(FloatMatrix * B)
{
	FloatMatrix Temp(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i, j) = this->at(i, j) - B->at(i, j);
		}
	}
	return Temp;
}


void FloatMatrix::Print()
{
	cout << "FloatMatrix " << setw(10) << m << "X" << n << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << setw(7) << this->at(i, j);
		}
		cout << endl;
	}
}


CSRMatrix::CSRMatrix()
{
	nRow = 0;
	nCol = 0;
	NonZero = 0;
	Values = NULL;
	RowIdx = NULL;
	ColIdx = NULL;
}


CSRMatrix::~CSRMatrix()
{
}


void CSRMatrix::Init(int NonZero, int nRow, int nCol, int* RowIdx, int*ColIdx)
{
	this->nRow = nRow;
	this->nCol = nCol;
	this->NonZero = NonZero;
	Values = new double[NonZero]();
	this->RowIdx = RowIdx;
	this->ColIdx = ColIdx;
}


// 返回i行j列的值
double & CSRMatrix::at(int i, int j)
{
	int idx;
	for (int iCol = RowIdx[i]; iCol < RowIdx[i + 1]; iCol++)
	{
		idx = ColIdx[iCol];
		if (idx == j)
		{
			return Values[iCol];
		}
	}
	double temp = 0;;
	return temp;
}


// 矩阵转置
CSRMatrix CSRMatrix::Trans(CSRMatrix * B)
{
	return CSRMatrix();
}


int CSRMatrix::GetRows()
{
	return nRow;
}


int CSRMatrix::GetCols()
{
	return nCol;
}


int CSRMatrix::GetNonZero()
{
	return NonZero;
}


// 打印数组
void CSRMatrix::Print()
{
	cout << "Compress Sparse Row Matrix:" << setw(10) << nRow << "X" << nCol << endl;
	for (int iRow = 0; iRow < nRow; iRow++)
	{
		for (int iCol = 0; iCol < nCol; iCol++)
		{
			cout << setw(7) << this->at(iRow, iCol);
		}
		cout << endl;
	}
	return void();
}
