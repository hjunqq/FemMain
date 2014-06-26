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
	// 构造大小为n的数组
	IntArray(int n);
	// 返回第i个函数值	
	int & at(int i);
	// 两个向量直接叠加
	IntArray Append(const IntArray & I)const;
	// 数组相加
	IntArray Add(const IntArray & I)const;
	// 数组复制
	IntArray Copy(const IntArray & I)const;
	IntArray(const IntArray & I);
	IntArray Cross(const IntArray &I)const;
	// 获取数组大小
	int GetSize();
	// 打印变量内容
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
	// 返回第i个元素
	double & at(int i);
	// 构造大小为n的数组
	FloatArray(int n);
	// 返回数组大小
	int GetSize();
	// 打印数组
	void Print();
	// 乘以标量
	FloatArray Times(double Scalar)const;
	// 点乘
	double Dot(FloatArray * );
	FloatArray Cross(FloatArray &F)const;
	// 清空
	int Clear();
	// 数组相加
	FloatArray Plus(const FloatArray & F)const;
	// 数组相减
	FloatArray Minus(const FloatArray & F )const;
	// 数组复制
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
	// 构造m行n列的矩阵
	FloatMatrix(int i, int j);
	FloatMatrix(const FloatMatrix & F);
	int GetSize();
	// 返回第m行n列的值
	double & at(int i, int j);
	// 转置
	FloatMatrix Trans() const;
	// 乘以数组
	FloatArray Mult(FloatArray & F )const;
	FloatMatrix Mult(double &D)const;
	// 乘以矩阵
	FloatMatrix Mult( FloatMatrix & F )const;
	// 矩阵相加
	FloatMatrix Plus(FloatMatrix & F)const;
	// 矩阵相减
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
	// 构造矩阵
	IntMatrix(int i, int j);
	IntMatrix(const IntMatrix &I);
	// 返回第m行n列元素
	int & at(int i, int j);
	double Determinant();
	IntMatrix  operator+(const IntMatrix & I)const;
	IntMatrix & operator=(const IntMatrix & I);
	IntMatrix  operator-(const IntMatrix & I)const;

};
 // Compressed Sparse Row Storage Format
 // 行压缩矩阵 
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
	// 返回i行j列的值
	double & at(int i, int j);
	// 矩阵转置
	CSRMatrix Trans(CSRMatrix & C)const;
	CSRMatrix Copy(CSRMatrix & C)const;
	int GetRows();
	int GetCols();
	int GetNonZero();
	// 打印数组
	void Print();
};
