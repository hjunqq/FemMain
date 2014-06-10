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
	// 构造大小为n的数组
	IntArray(int n);
	// 返回第i个函数值	
	int & at(int i);
	// 两个向量直接叠加
	IntArray Append( IntArray & TArray);
	// 数组相加
	IntArray Add(IntArray & TArrray);
	// 数组复制
	IntArray Copy(IntArray *TArray);
	// 获取数组大小
	int GetSize();
	// 打印变量内容
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
	// 返回第i个元素
	double & at(int i);
	// 构造大小为n的数组
	FloatArray(int n);
	// 返回数组大小
	int GetSize();
	// 打印数组
	void Print();
	// 乘以标量
	FloatArray Times(double Scalar);
	// 点乘
	double Dot(FloatArray * );
	// 清空
	int Clear();
	// 数组相加
	FloatArray Plus(FloatArray * );
	// 数组相减
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
	// 构造m行n列的矩阵
	FloatMatrix(int i, int j);
	// 返回第m行n列的值
	double & at(int i, int j);
	// 转置
	FloatMatrix Trans();
	// 乘以数组
	FloatArray Mult(FloatArray * );
	// 乘以矩阵
	FloatMatrix Mult(FloatMatrix * );
	// 矩阵相加
	FloatMatrix Plus(FloatMatrix *);
	// 矩阵相减
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
	// 构造矩阵
	IntMatrix(int i, int j);
	// 返回第m行n列元素
	int & at(int i, int j);
};
 // Compressed Sparse Row Storage Format
 // 行压缩矩阵 
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
	// 返回i行j列的值
	double & at(int i, int j);
	// 矩阵转置
	CSRMatrix Trans(CSRMatrix * B);
	int GetRows();
	int GetCols();
	int GetNonZero();
	// 打印数组
	void Print();
};
