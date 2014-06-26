#pragma once
#include <fstream> 
#include <sstream>
#include <iomanip>
#include <iostream>
using namespace std;
class Time
{
private:
	int hours;
	int minutes;
public:
	Time();
	Time(int h, int m = 0);
	void AddMin(int m);
	void AddHr(int h);
	void Reset(int h = 0, int m = 0);
	Time operator+(const Time &t)const;
	Time& operator=(const Time &t);
	void Show() const;
};
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
	IntArray Append(const IntArray & I)const;
	// �������
	IntArray Add(const IntArray & I)const;
	// ���鸴��
	IntArray Copy(const IntArray & I)const;
	IntArray(const IntArray & I);
	IntArray Cross(const IntArray &I)const;
	// ��ȡ�����С
	int GetSize();
	// ��ӡ��������
	void Print();
	void Set(int Value);
	void SetSize(int Size);
	

	IntArray operator+(const IntArray & I)const;
	IntArray & operator=(const IntArray & I);
	IntArray operator-(const IntArray & I)const;
};

class FloatArray
{
public:
	FloatArray();
	virtual ~FloatArray();
	bool IsNULL();
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
	FloatArray Times(double Scalar)const;
	// ���
	double Dot(FloatArray * );
	FloatArray Cross(FloatArray &F)const;
	// ���
	int Clear();
	// �������
	FloatArray Plus(const FloatArray & F)const;
	// �������
	FloatArray Minus(const FloatArray & F )const;
	// ���鸴��
	FloatArray Copy(const FloatArray & F)const;
	FloatArray(const FloatArray & F);
	void SetSize(int Size);
	FloatArray  operator+(const FloatArray & I)const;
	FloatArray & operator=(const FloatArray & I);
	FloatArray  operator-(const FloatArray & I)const;
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
	int GetSize();
	// ���ص�m��n�е�ֵ
	double & at(int i, int j);
	// ת��
	FloatMatrix Trans() const;
	// ��������
	FloatArray Mult(FloatArray & F )const;
	FloatMatrix Mult(double &D)const;
	// ���Ծ���
	FloatMatrix Mult( FloatMatrix & F )const;
	// �������
	FloatMatrix Plus(FloatMatrix & F)const;
	// �������
	FloatMatrix Minus(FloatMatrix & F )const;
	FloatMatrix Copy(FloatMatrix & F)const;
	void SetSize(int m, int n);
	double Determinant();

	FloatMatrix Inverse();

	void Print();

	FloatMatrix  operator+(const FloatMatrix & I)const;
	FloatMatrix  & operator=(const FloatMatrix & I);
	FloatMatrix  operator-(const FloatMatrix & I)const;
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
	double Determinant();
	IntMatrix  operator+(const IntMatrix & I)const;
	IntMatrix & operator=(const IntMatrix & I);
	IntMatrix  operator-(const IntMatrix & I)const;

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
	CSRMatrix Trans(CSRMatrix & C)const;
	CSRMatrix Copy(CSRMatrix & C)const;
	int GetRows();
	int GetCols();
	int GetNonZero();
	// ��ӡ����
	void Print();
};
