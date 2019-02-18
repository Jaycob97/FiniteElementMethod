#pragma once
#include "Element.h"

class Jacobian
{
public:
	double N[4][4], ksi[4], eta[4], InterpolatedCoordinates[4][2], dNdksi[4][4], dNdeta[4][4], jacobian[4][4], detJacobian[4], inversedJacobian[4][4];

	Jacobian();
	~Jacobian();
	void calculateShapeFunctions();
	void calculateInterpolatedCoordinates(Element element);
	void calculateShapeFunctionsDerivatives();
	void calculateJacobian(Element element);
};