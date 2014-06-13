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
	IntArray Append(const IntArray & I);
	// �������
	IntArray Add(const IntArray & I);
	// ���鸴��
	IntArray Copy(const IntArray & I);
	IntArray(const IntArray & I);
	// ��ȡ�����С
	int GetSize();
	// ��ӡ��������
	void Print();
	void Set(int Value);

	IntArray & operator+(const IntArray & I);
	IntArray & operator=(const IntArray & I);
	IntArray & operator-(const IntArray & I);
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
	FloatArray Plus(const FloatArray & F);
	// �������
	FloatArray Minus(const FloatArray & F );
	// ���鸴��
	FloatArray Copy(const FloatArray & F);
	FloatArray(const FloatArray & F);

	FloatArray & operator+(const FloatArray & I);
	FloatArray & operator=(const FloatArray & I);
	FloatArray & operator-(const FloatArray & I);
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
	FloatMatrix(const FloatMatrix & F);
	// ���ص�m��n�е�ֵ
	double & at(int i, int j);
	// ת��
	FloatMatrix Trans();
	// ��������
	FloatArray Mult(FloatArray & F );
	FloatMatrix Mult(double &D);
	// ���Ծ���
	FloatMatrix Mult( FloatMatrix & F );
	// �������
	FloatMatrix Plus(FloatMatrix & F);
	// �������
	FloatMatrix Minus(FloatMatrix & F );
	FloatMatrix Copy(FloatMatrix & F);
	
	double Determinant();

	FloatMatrix Inverse();

	void Print();

	FloatMatrix & operator+(const FloatMatrix & I);
	FloatMatrix & operator=(const FloatMatrix & I);
	FloatMatrix & operator-(const FloatMatrix & I);
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
	IntMatrix(const IntMatrix &I);
	// ���ص�m��n��Ԫ��
	int & at(int i, int j);

	IntMatrix & operator+(const IntMatrix & I);
	IntMatrix & operator=(const IntMatrix & I);
	IntMatrix & operator-(const IntMatrix & I);

};
 // Compressed Sparse Row Storage Format
 // ��ѹ������ 
class CSRMatrix
{
public:
	CSRMatrix();
	CSRMatrix(CSRMatrix &C);
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
	CSRMatrix Trans(CSRMatrix & C);
	CSRMatrix Copy(CSRMatrix & C);
	int GetRows();
	int GetCols();
	int GetNonZero();
	// ��ӡ����
	void Print();
};
