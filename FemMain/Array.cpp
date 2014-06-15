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

void IntArray::Set(int Value)
{
	for (int i = 0; i < size; i++)
	{
		Values[i] = Value;
	}
}

// 两个向量直接叠加
IntArray IntArray::Append(const IntArray & I)
{
	IntArray Temp;
	int nNew;
	nNew = size + I.size;
	Temp.Values = new int[nNew]();
	Temp.size = nNew;
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = Values[i];
	}
	for (int i = size; i < nNew; i++)
	{
		Temp.Values[i] = I.Values[i - size];
	}
	return Temp;
}


// 数组相加
IntArray IntArray::Add(const IntArray & I)
{
	IntArray Temp;
	Temp.Values = new int[size]();
	Temp.size =size;
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = Values[i] + I.Values[i];
	}
	return Temp;
}

// 数组复制
IntArray::IntArray(const IntArray & I)
{
	size = I.size;
	Values = new int[size];
	for (int i = 0; i < size; i++)
	{
		Values[i] = I.Values[i];
	}
}
IntArray IntArray::Copy(const IntArray & I)
{
	IntArray Temp;
	Temp.size = I.size;
	Temp.Values = new int[size];
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = I.Values[i];
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

IntArray &  IntArray::operator + (const IntArray & I)
{
	IntArray T;
	T.size = I.size;
	T.Values = new int[T.size];
	for (int i = 0; i < T.size; i++)
	{
		T.Values[i] =Values[i]+ I.Values[i];
	}
	return T;
}
IntArray & IntArray::operator = (const IntArray & I)
{
	IntArray T;
	T.size =  I.size;
	T.Values = new int[T.size];
	for (int i = 0; i < T.size; i++)
	{
		T.Values[i] = I.Values[i];
	}
	return T;
}
IntArray & IntArray::operator - (const IntArray & I)
{
	IntArray T;
	T.size = I.size;
	T.Values = new int[T.size];
	for (int i = 0; i < T.size; i++)
	{
		T.Values[i] =Values[i]- I.Values[i];
	}
	return T;
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
	for (int i = 0; i < size; i++)
	{
		cout << setw(15) << Values[i] ;
	}
	cout << endl;
}


// 乘以标量
FloatArray FloatArray::Times(double Scalar)
{
	FloatArray TArray;
	TArray.size = size;
	TArray.Values = new double[TArray.size];
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
	for (int i = 0; i < size; i++)
	{
		Values[i] = 0;
	}
	return 0;
}


// 数组相加
FloatArray FloatArray::Plus(const FloatArray & F)
{
	FloatArray Temp;
	Temp.size = size;
	Temp.Values = new double[Temp.size];
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = Values[i] + F.Values[i];
	}
	return Temp;
}


// 数组相减
FloatArray FloatArray::Minus(const FloatArray & F)
{
	FloatArray Temp;
	Temp.size = size;
	Temp.Values = new double[Temp.size];
	for (int i = 0; i < size; i++)
	{
		Temp.Values[i] = Values[i] - F.Values[i];
	}
	return Temp;
}

 // 数组复制
FloatArray FloatArray::Copy(const FloatArray & F)
{
	FloatArray Temp;
	Temp.size = F.size;
	Temp.Values = new double[Temp.size]();
	for (int i = 0; i < Temp.size; i++)
	{
		Temp.Values[i] = F.Values[i];
	}
	return Temp;
}
FloatArray::FloatArray(const FloatArray & F)
{
	size = F.size;
	Values = new double[size]();
	for (int i = 0; i < size; i++)
	{
		Values[i] = F.Values[i];
	}
}

FloatArray & FloatArray::operator + (const FloatArray & I)
{
	FloatArray T;
	T.size = I.size;
	T.Values = new double[T.size];
	for (int i = 0; i < T.size; i++)
	{
		T.Values[i] =this->Values[i]+ I.Values[i];
	}
	return T;
}
FloatArray & FloatArray::operator = (const FloatArray & I)
{
	FloatArray T;
	T.size = I.size;
	T.Values = new double[T.size];
	for (int i = 0; i < T.size; i++)
	{
		T.Values[i] = I.Values[i];
	}
	return T;
}
FloatArray & FloatArray::operator - (const FloatArray & I)
{
	FloatArray T;
	return T;
}

FloatMatrix::FloatMatrix()
{
	m = 0;
	n = 0;
	Values = new double[m*n]();
}

FloatMatrix::FloatMatrix(const FloatMatrix &F)
{
	m = F.m;
	n = F.n;
	Values = new double[m*n]();
	for (int i = 0; i < m*n; i++)
	{
		Values[i] = F.Values[i];
	}
}

FloatMatrix::~FloatMatrix()
{
}


IntMatrix::IntMatrix()
{
	m = 0;
	n = 0;
	Values = new int[m*n]();
}

IntMatrix::IntMatrix(const IntMatrix & I)
{
	m = I.m;
	n = I.n;
	Values = new int[m*n]();
	for (int i = 0; i < m*n; i++)
	{
		Values[i] = I.Values[i];
	}
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
IntMatrix & IntMatrix::operator + (const IntMatrix & I)
{
	IntMatrix T;
	T.m = m;
	T.n = n;
	T.Values = new int[m + n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			T.at(i, j) = this->at(i, j) + I.Values[i*n+j];
		}
	}
	return T;
}
IntMatrix & IntMatrix::operator=(const IntMatrix & I)
{
	IntMatrix T;
	T.m = m;
	T.n = n;
	T.Values = new int[m + n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			T.at(i, j) = I.Values[i*n + j];
		}
	}
	return T;
}
IntMatrix & IntMatrix::operator-(const IntMatrix & I)
{
	IntMatrix T;
	T.m = m;
	T.n = n;
	T.Values = new int[m + n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			T.at(i, j) = this->at(i, j) - I.Values[i*n + j];
		}
	}
	return T;
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

FloatMatrix FloatMatrix::Mult(double &D)
{
	FloatMatrix Temp(m,n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i,j) = this->at(i, j)*D;
		}
	}
	return Temp;
}
// 乘以数组
FloatArray FloatMatrix::Mult(FloatArray & F)
{
	FloatArray Temp(m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i) += F.at(i) * this->at(i, j);
		}
	}
	return Temp;
}

 // 乘以矩阵
FloatMatrix FloatMatrix::Mult( FloatMatrix & F)
{
	FloatMatrix Temp(m, F.n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < F.n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				Temp.at(i, j) += this->at(i, k)*F.at(k,j);
			}			
		}
	}
	return Temp;
}

 // 矩阵相加
FloatMatrix FloatMatrix::Plus(FloatMatrix & F)
{
	FloatMatrix Temp(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i, j) = this->at(i, j) + F.at(i, j);
		}
	}
	return Temp;
}

 // 矩阵相减
FloatMatrix FloatMatrix::Minus(FloatMatrix & F)
{
	FloatMatrix Temp(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Temp.at(i, j) = this->at(i, j) - F.at(i, j);
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
			cout << setw(15) << this->at(i, j);
		}
		cout << endl;
	}
}

double FloatMatrix::Determinant()
{
	double Det=0;
	if (this->m == 1)
	{
		return this->at(0, 0);
	}
	else
	{
		FloatMatrix *T;
		T = new FloatMatrix(this->m - 1, this->n - 1);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				for (int k = 0; k < m; k++)
				{
					if (j < i && k < i)
					{
						T->at(j, k) = this->at(j, k);
					}
					if (j>i && k < i)
					{
						T->at(j - 1, k) = this->at(j, k);
					}
					if (j<i && k>i)
					{
						T->at(j, k - 1) = this->at(j, k);
					}
					if (j>i && k > i)
					{
						T->at(j - 1, k - 1) = this->at(j, k);
					}
				}
			}

			Det += this->at(i, 0)* pow(-1, m)*T->Determinant();
		}
		return Det;
	}
}

FloatMatrix FloatMatrix::Inverse()
{
	FloatMatrix T1(m,n),T2(m-1,n-1);
	double det;
	det = this->Determinant();
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < m; k++)
			{
				for (int l = 0; l < m; l++)
				{
					if (k < i && l < j)
					{
						T2.at(k, l) = this->at(k, l);
					}
					if (k>i && l < j)
					{
						T2.at(k - 1, l) = this->at(k, l);
					}
					if (k<i && l>j)
					{
						T2.at(k, l - 1) = this->at(k, l);
					}
					if (k>i && l > j)
					{
						T2.at(k - 1, l - 1) = this->at(k, l);
					}
				}
			}
			T1.at(i, j) = T2.Determinant()/det;
		}
	}
	return T1;
}


FloatMatrix & FloatMatrix::operator + (const FloatMatrix & I)
{
	FloatMatrix T;
	T.m = m;
	T.n = n;
	T.Values = new double[m + n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			T.at(i, j) = this->at(i, j) + I.Values[i*n + j];
		}
	}
	return T;
}
FloatMatrix & FloatMatrix::operator = (const FloatMatrix & I)
{
	FloatMatrix T;
	T.m = m;
	T.n = n;
	T.Values = new double[m + n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			T.at(i, j) =I.Values[i*n + j];
		}
	}
	return T;
}
FloatMatrix &  FloatMatrix::operator - (const FloatMatrix & I)
{
	FloatMatrix T;
	T.m = m;
	T.n = n;
	T.Values = new double[m + n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			T.at(i, j) = this->at(i, j) - I.Values[i*n + j];
		}
	}
	return T;
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

CSRMatrix::CSRMatrix(CSRMatrix &C)
{
	nRow = C.nRow;
	nCol = C.nCol;
	NonZero = C.NonZero;
	Values = new double[NonZero]();
	RowIdx = new int [nRow + 1]();
	ColIdx = new int[NonZero]();
	for (int i = 0; i < NonZero; i++)
	{
		ColIdx[i] = C.ColIdx[i];
		Values[i] = C.Values[i];
	}
	for (int i = 0; i < nRow+1; i++)
	{
		RowIdx[i] = C.RowIdx[i];
	}
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
	Zero = 0;
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
	return Zero;
}


// 矩阵转置
CSRMatrix CSRMatrix::Trans(CSRMatrix &C)
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
