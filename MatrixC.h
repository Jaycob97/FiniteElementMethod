#pragma once
#include "Jacobian.h"

class MatrixC
{
public:
    double C[4][4];

    MatrixC();
    ~MatrixC();

    void calculateMatrixC(Jacobian jacobian, double c, double ro);
};