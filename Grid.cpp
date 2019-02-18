#ifdef _MSC_VER
#endif

#include "stdafx.h"
#include "Grid.h"
#include "MatrixH.h"
#include "MatrixC.h"
#include "MatrixHBC.h"
#include "VectorP.h"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <limits>

Grid::Grid()
{
}


Grid::~Grid()
{
}

void Grid::generateGrid(std::string file)
{
	std::ifstream stream(file);
	if (stream)
	{
	stream >> H;
	stream >> L;
	stream >> nH;
	stream >> nL;
	stream >> k;
	stream >> c;
	stream >> ro;
	stream >> alpha;
	stream >> ambientTemperature;
	stream >> simulationStepTime;
	stream >> initialTemperature;
	stream >> simulationTime;
	}
	else
	{
		std::cout << "Blad odczytu" << std::endl;
		exit(0);
	}
    std::cout << "H: " << H << std::endl;
    std::cout << "L: " << L << std::endl;
    std::cout << "nH: " << nH << std::endl;
    std::cout << "nL: " << nL << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "ro: " << ro << std::endl;
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "ambientTemperature: " << ambientTemperature << std::endl;
    std::cout << "simulationStepTime: " << simulationStepTime << std::endl;
    std::cout << "initialTemperature: " << initialTemperature << std::endl;
    std::cout << "simulationTime: " << simulationTime << std::endl;

    nodeList = new Node[nH*nL];

    for (int j = 0; j < nL; j++)
    {
        for (int i = 0; i < nH; i++)
        {
            nodeList[i + (j*nH)].y = i * (H / (nH - 1));
            nodeList[i + (j*nH)].x = j * (L / (nL - 1));
            nodeList[i + (j*nH)].t = initialTemperature;
        }
    }

    elementList = new Element[(nH - 1)*(nL - 1)];

    for (int j = 0; j < (nL - 1); j++)
    {
        for (int i = 0; i < (nH - 1); i++)
        {
            elementList[i + (j*(nH - 1))].nodeID[0] = nodeList[i + (j * nH)];
            elementList[i + (j*(nH - 1))].nodeID[1] = nodeList[i + (j * nH) + nH];
            elementList[i + (j*(nH - 1))].nodeID[2] = nodeList[i + (j * nH) + nH + 1];
            elementList[i + (j*(nH - 1))].nodeID[3] = nodeList[i + (j*nH) + 1];

            elementList[i + (j*(nH - 1))].ID[0] = (i + (j * nH));
            elementList[i + (j*(nH - 1))].ID[1] = (i + (j * nH) + nH);
            elementList[i + (j*(nH - 1))].ID[2] = (i + (j * nH) + nH + 1);
            elementList[i + (j*(nH - 1))].ID[3] = (i + (j*nH) + 1);
        }
    }

    // warunki brzegowe
    for(int i = 0; i < nH-1; i++)
    {
        elementList[i].isSurfaceHeated[3] = true;
    }

    for(int i = ((nH-1)*(nL-1))-1; i > ((nH-1)*(nL-1))-1-(nH-1); i--)
    {
        elementList[i].isSurfaceHeated[1] = true;
    }

    for(int i = 0; i < nL-1; i++)
    {
        elementList[i*(nH-1)].isSurfaceHeated[0] = true;
    }

    for(int i = 0; i < nL-1; i++)
    {
        elementList[i*(nH-1)+(nH-2)].isSurfaceHeated[2] = true;
    }

    fclose(stdin);


}

void Grid::aggregateGrid() {
    double *temp;

	temp = new double[nL*nH];

    matrixHAndVectorP = new double*[nL*nH];
    for(int i = 0; i < nL*nH; ++i)
        matrixHAndVectorP[i] = new double[nL*nH];

    globalMatrixH = new double*[nL*nH];
    for(int i = 0; i < nL*nH; ++i)
        globalMatrixH[i] = new double[nL*nH];

    globalMatrixHBC = new double*[nL*nH];
    for(int i = 0; i < nL*nH; ++i)
        globalMatrixHBC[i] = new double[nL*nH];

    globalMatrixC = new double*[nL*nH];
    for(int i = 0; i < nL*nH; ++i)
        globalMatrixC[i] = new double[nL*nH];

    globalVectorP = new double[nL*nH];

    vectorP = new double[nL*nH];

    for(int i = 0; i < nL*nH; i++)
    {
        globalVectorP[i]=0;

        for(int j = 0; j < nH*nH; j++)
        {
            globalMatrixH[i][j] = 0;
            globalMatrixHBC[i][j] = 0;
            globalMatrixC[i][j] = 0;
        }
    }

    for(int i = 0; i < (nH-1)*(nL-1); i++) {
        Jacobian tempJacobian;

        tempJacobian.calculateInterpolatedCoordinates(elementList[i]);
        tempJacobian.calculateShapeFunctionsDerivatives();
        tempJacobian.calculateJacobian(elementList[i]);

        MatrixH tempMatrixH;
        tempMatrixH.calculateMatrixH(tempJacobian, k);
        for(int k=0; k<4; k++)
        {
            for(int j = 0; j<4; j++)
            {
                globalMatrixH[elementList[i].ID[k]][elementList[i].ID[j]] += tempMatrixH.H[k][j];
            }
        }

        MatrixC tempMatrixC;
        tempMatrixC.calculateMatrixC(tempJacobian, c, ro);
        for(int k=0; k<4; k++)
        {
            for(int j = 0; j<4; j++)
            {
                globalMatrixC[elementList[i].ID[k]][elementList[i].ID[j]] += tempMatrixC.C[k][j];
            }
        }

        MatrixHBC tempMatrixHBC;
        tempMatrixHBC.alpha = alpha;
        tempMatrixHBC.CalculateMatrixHBC(elementList[i]);
        for(int k=0; k<4; k++)
        {
            for(int j = 0; j<4; j++)
            {
                globalMatrixHBC[elementList[i].ID[k]][elementList[i].ID[j]] += tempMatrixHBC.matrixH[k][j];
            }
        }

        vec = new double[nL*nH];
        for(int i=0; i<nL*nH; i++)
        {
            vec[i]=initialTemperature;
        }

        for(int i=0; i<nL*nH; i++)
        {
            temp[i]=0;
        }

        for(int i = 0; i < nL*nH; i++){
            for(int j = 0; j < nL*nH; j++)
            {
                temp[i] += ((globalMatrixC[i][j]/simulationStepTime) * vec[j]);
            }
        }

        VectorP tempVectorP;
        tempVectorP.alpha = alpha;
        tempVectorP.ambientTemperature = ambientTemperature;
        tempVectorP.CalculateVectorP(elementList[i]);
        for(int k=0; k<4; k++)
        {
                globalVectorP[elementList[i].ID[k]] += tempVectorP.vectorP[k];
        }
    }
    for (int j = 0; j < nL*nH; j++) {

        vectorP[j] = globalVectorP[j];
        globalVectorP[j] = -globalVectorP[j] + temp[j];
        for (int k = 0; k < nL*nH; k++) {
            globalMatrixH[j][k] += globalMatrixHBC[j][k]+ (globalMatrixC[j][k]/simulationStepTime);

        }
    }
}

void Grid::setValueOfMatrixHAndVectorP() {
    for (int i = 0; i < (nH * nL); ++i) {
        for (int j = 0; j < (nH * nL) + 1; ++j) {
            if (j == (nH * nL)) {
                matrixHAndVectorP[i][j] = globalVectorP[i];
                continue;
            }
            matrixHAndVectorP[i][j] = globalMatrixH[i][j];
        }
    }
}

void Grid::calculateTemperatures() {
    int n = static_cast<int>(nH * nL);
    for (int j = 0; j < simulationTime; j = j + simulationStepTime) {
        double max = std::numeric_limits<double>::min();
        double min = std::numeric_limits<double>::max();
        setValueOfMatrixHAndVectorP();
        gaussMethod(n);
        for(int i = 0; i < (nL*nH); i++)
        {
            if(vec[i] > max)
                max = vec[i];
            if(vec[i] < min)
                min = vec[i];
        }
        std::cout << j+simulationStepTime <<"s"<<" "<< "minTemp " << std::fixed << min << " " << "maxTemp " << " " << max << std::endl;
        updateVectorP();
    }
}

void Grid::updateVectorP() {
    double val1 = 0;
    for (int j = 0; j < (nL * nH); ++j) {
        globalVectorP[j] = 0;
        for (int i = 0; i < (nL * nH); ++i) {
            val1 += ((globalMatrixC[j][i] / simulationStepTime) * vec[i]);
        }
        globalVectorP[j] = -vectorP[j] + val1;
        val1 = 0;
    }
}

bool Grid::gaussMethod(int n) {
    const double eps = 1e-12;
    int i, j, k;
    double m, s;
    for(int i = 0; i < (nL*nH); i++)
    {
        vec[i] = 0;
    }

    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (fabs(matrixHAndVectorP[i][i]) < eps) return false;
            m = -matrixHAndVectorP[j][i] / matrixHAndVectorP[i][i];
            for (k = i + 1; k <= n; k++)
                matrixHAndVectorP[j][k] += m * matrixHAndVectorP[i][k];
        }
    }

    for (i = n - 1; i >= 0; i--) {
        s = matrixHAndVectorP[i][n];
        for (j = n - 1; j >= i + 1; j--)
            s -= matrixHAndVectorP[i][j] * vec[j];
        if (fabs(matrixHAndVectorP[i][i]) < eps) return false;
        vec[i] = s / matrixHAndVectorP[i][i];
    }
    return true;
}
