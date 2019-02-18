#include "stdafx.h"
#include "Jacobian.h"
#include <iostream>
#include <cmath>

Jacobian::Jacobian()
{
	ksi[0] = -1 / sqrt(3);
	ksi[1] = -ksi[0];
	ksi[2] = ksi[1];
	ksi[3] = -ksi[2];

	eta[0] = ksi[0];
	eta[1] = eta[0];
	eta[2] = -eta[1];
	eta[3] = eta[2];

	calculateShapeFunctions();
}

Jacobian::~Jacobian()
{
}

void Jacobian::calculateShapeFunctions()
{
	for(int i = 0; i < 4; i++)
	{
		N[i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		N[i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		N[i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		N[i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}
}

void Jacobian::calculateInterpolatedCoordinates(Element element)
{
	for(int i = 0; i < 4; i++)
	{
		InterpolatedCoordinates[i][0] = N[i][0] * element.nodeID[0].x + N[i][1] * element.nodeID[1].x + N[i][2] * element.nodeID[2].x + N[i][3] * element.nodeID[3].x;
		InterpolatedCoordinates[i][1] = N[i][0] * element.nodeID[0].y + N[i][1] * element.nodeID[1].y + N[i][2] * element.nodeID[2].y + N[i][3] * element.nodeID[3].y;
	}
}

void Jacobian::calculateShapeFunctionsDerivatives()
{
	for (int i = 0; i < 4; i++)
	{
		dNdksi[0][i] = -0.25*(1 - eta[i]);
		dNdksi[1][i] = 0.25*(1 - eta[i]);
		dNdksi[2][i] = 0.25*(1 + eta[i]);
		dNdksi[3][i] = -0.25*(1 + eta[i]);

		dNdeta[0][i] = -0.25*(1 - ksi[i]);
		dNdeta[1][i] = -0.25*(1 + ksi[i]);
		dNdeta[2][i] = 0.25*(1 + ksi[i]);
		dNdeta[3][i] = 0.25*(1 - ksi[i]);
	}

}

void Jacobian::calculateJacobian(Element element)
{
	for(int i = 0; i < 4; i++)
	{
		jacobian[i][0] = element.nodeID[0].x * dNdksi[0][i] + element.nodeID[1].x * dNdksi[1][i] + element.nodeID[2].x * dNdksi[2][i] + element.nodeID[3].x * dNdksi[3][i];
		jacobian[i][1] = element.nodeID[0].y * dNdksi[0][i] + element.nodeID[1].y * dNdksi[1][i] + element.nodeID[2].y * dNdksi[2][i] + element.nodeID[3].y * dNdksi[3][i];
		jacobian[i][2] = element.nodeID[0].x * dNdeta[0][i] + element.nodeID[1].x * dNdeta[1][i] + element.nodeID[2].x * dNdeta[2][i] + element.nodeID[3].x * dNdeta[3][i];
		jacobian[i][3] = element.nodeID[0].y * dNdeta[0][i] + element.nodeID[1].y * dNdeta[1][i] + element.nodeID[2].y * dNdeta[2][i] + element.nodeID[3].y * dNdeta[3][i];

		detJacobian[i] = jacobian[i][0] * jacobian[i][3] - jacobian[i][1] * jacobian[i][2];

		inversedJacobian[i][0] = jacobian[i][3] / detJacobian[i];
		inversedJacobian[i][1] = -(jacobian[i][1] / detJacobian[i]);
		inversedJacobian[i][2] = -(jacobian[i][2] / detJacobian[i]);
		inversedJacobian[i][3] = jacobian[i][0] / detJacobian[i];
	}
}
