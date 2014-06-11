#include "Array.h"

IntArray::IntArray()
{
	size = 0;
}


IntArray::~IntArray()
{
}


// �����СΪn������
IntArray::IntArray(int n)
{
	size = n;
	Values = new int[n]();
}


 //���ص�i������ֵ	
int & IntArray::at(int i)
{
	return Values[i];
}



// ��������ֱ�ӵ���
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


// �������
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

// ���鸴��
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
// ��ȡ�����С
int IntArray::GetSize()
{
	return size;
}


// ��ӡ��������
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


// ���ص�i��Ԫ��
double & FloatArray::at(int i)
{
	return Values[i];
}


// �����СΪn������
FloatArray::FloatArray(int n)
{
	size = n;
	Values = new double[n]();
}


// ���������С
int FloatArray::GetSize()
{
	return size;
}


// ��ӡ����
void FloatArray::Print()
{
	for (int i = 0; i < size; i++)
	{
		cout << setw(15) << Values[i] ;
	}
	cout << endl;
}


// ���Ա���
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


// ���
double FloatArray::Dot(FloatArray * B)
{
	double sum=0.0;
	for (int i = 0; i < size; i++)
	{
		sum += Values[i] * B->at(i);
	}
	return sum;
}


// ���
int FloatArray::Clear()
{
	delete[] Values;
	size = 0;
	return 0;
}


// �������
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


// �������
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

 // ���鸴��
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


// �������
IntMatrix::IntMatrix(int i, int j)
{
	m = i; n = j;
	Values = new int[m*n]();
}


// ���ص�m��n��Ԫ��
int & IntMatrix::at(int i, int j)
{
	return Values[i*n + j];
}


// ����m��n�еľ���
FloatMatrix::FloatMatrix(int i, int j)
{
	m = i; n = j;
	Values = new double[m*n]();
}


// ���ص�m��n�е�ֵ
double & FloatMatrix::at(int i, int j)
{
	return Values[i*n + j];
}


// ת��
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


// ��������
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

 // ���Ծ���
FloatMatrix FloatMatrix::Mult( FloatMatrix & F)
{
	FloatMatrix Temp(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < F.n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				Temp.at(i, j) = this->at(i, k)*F.at(k,j);
			}			
		}
	}
	return Temp;
}

 // �������
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

 // �������
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


// ����i��j�е�ֵ
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


// ����ת��
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


// ��ӡ����
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
