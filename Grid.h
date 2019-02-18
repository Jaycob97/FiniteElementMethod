#pragma once
#include "Element.h"
#include "MatrixH.h"
#include "MatrixHBC.h"
#include "MatrixC.h"
#include "VectorP.h"
#include <string>
#include <fstream>

class Grid
{
public:
	Node * nodeList;
	Element* elementList;
	double **globalMatrixH;
	double **globalMatrixHBC;
	double **globalMatrixC;
	double *globalVectorP;
	float H, L;
	int nH, nL;
	double k, c, ro, alpha, ambientTemperature, simulationTime, simulationStepTime, initialTemperature;
	double *vec;
	double **matrixHAndVectorP;
	double *vectorP;

	Grid();
	~Grid();
	void generateGrid(std::string file);
	void aggregateGrid();
	void updateVectorP();
	void calculateTemperatures();
	bool gaussMethod(int n);
	void setValueOfMatrixHAndVectorP();
};
