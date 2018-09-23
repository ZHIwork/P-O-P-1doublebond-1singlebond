#pragma once
#include <iostream>
#include "CONSTANT.h"
#include "EnergyForce.h"
#include "Tool.h"
#include "mkl.h"

class BFGS
{
private:
	//all of these arrays are in Cartesian Coordinates
	int dimension;
	T inverseB0[193 * 3 * 193 * 3];
	T x_before[193*3], x_after[193 * 3];
	T inverseB_before[193 * 3 * 193 * 3], inverseB_after[193 * 3 * 193 * 3];
	T g_before[193 * 3], g_after[193 * 3];
	T S_before[193 * 3];
	T theta[193 * 3], y[193 * 3];
	T Alpha, Beta, Gamma, t;
	T temp2Dimension[193 * 3 * 193 * 3], tempDimension[193 * 3], tempOne[1], I[193 * 3 * 193 * 3];
	T XCart[193 * 3];
	T XDirect[193 * 3];
	T gX[193 * 3];
public:
	BFGS();
	~BFGS();

public:
	void LineSearch(T &totalEnergy, const int doubleO[2][2]);
	T GetAlpha(const T &alpha, const T Sk[], const T gk[], const T &totalEnergy, const int doubleO[2][2]);
	void setZero(T a[], const int &N);
	void assignValuefor1Darray(const T a[], T b[], const int &N);
	void add(const T a[], const T b[], T c[], const int &N);
	void subtract(const T a[], const T b[], T c[], const int &N);
	void numMultiply(const T &num, const T a[], T b[], const int &N);
	T specialMultiply(const T a[], const T b[], const int &N);
	void convertDirecttoCartesianCoordinates(const T coorDirect[], T coorCartesian[], const int &N);
	void convertandLimitCarttoDirectCoordinates(const T coorCartesian[], T coorDirect[], const int &N);
	void limitCoor(T a[], const int &N);
};

