#pragma once
#include <fstream> 
#include <sstream>
#include <iomanip>
#include <iostream>


using namespace std;

class IntArray
{
public:
	IntArray();
	virtual ~IntArray();
private:
	int size;
	int *Values;
public:
	// �����СΪn������
	IntArray(int n);
	// ���ص�i������ֵ	
	int & at(int i);
	// ��������ֱ�ӵ���
	IntArray Append( IntArray & TArray);
	// �������
	IntArray Add(IntArray & TArrray);
	// ���鸴��
	IntArray Copy(IntArray *TArray);
	// ��ȡ�����С
	int GetSize();
	// ��ӡ��������
	void Print();
};

class FloatArray
{
public:
	FloatArray();
	virtual ~FloatArray();
private:
	int size;
	double *Values;
public:
	// ���ص�i��Ԫ��
	double & at(int i);
	// �����СΪn������
	FloatArray(int n);
	// ���������С
	int GetSize();
	// ��ӡ����
	void Print();
	// ���Ա���
	FloatArray Times(double Scalar);
	// ���
	double Dot(FloatArray * );
	// ���
	int Clear();
	// �������
	FloatArray Plus(FloatArray * );
	// �������
	FloatArray Minus(FloatArray * );
};

class FloatMatrix
{
public:
	FloatMatrix();
	virtual ~FloatMatrix();
private:
	int m;
	int n;
	double *Values;
public:
	// ����m��n�еľ���
	FloatMatrix(int i, int j);
	// ���ص�m��n�е�ֵ
	double & at(int i, int j);
	// ת��
	FloatMatrix Trans();
	// ��������
	FloatArray Mult(FloatArray * );
	// ���Ծ���
	FloatMatrix Mult(FloatMatrix * );
	// �������
	FloatMatrix Plus(FloatMatrix *);
	// �������
	FloatMatrix Minus(FloatMatrix * );
	void Print();
};

class IntMatrix
{
public:
	IntMatrix();
	virtual ~IntMatrix();
private:
	int m;
	int n;
	int * Values;
public:
	// �������
	IntMatrix(int i, int j);
	// ���ص�m��n��Ԫ��
	int & at(int i, int j);
};
 // Compressed Sparse Row Storage Format
 // ��ѹ������ 
class CSRMatrix
{
public:
	CSRMatrix();
	~CSRMatrix();
private:
	int nRow;
	int nCol;
	int NonZero;
	double *Values;
	int *RowIdx;
	int *ColIdx;
	double Zero;
public:
	void Init(int NonZero, int nRow, int nCol, int* RowIdx,int*ColIdx);
	// ����i��j�е�ֵ
	double & at(int i, int j);
	// ����ת��
	CSRMatrix Trans(CSRMatrix * B);
	int GetRows();
	int GetCols();
	int GetNonZero();
	// ��ӡ����
	void Print();
};
