#ifndef FINITEELEMENTMETHOD_MATRIXHBC_H
#define FINITEELEMENTMETHOD_MATRIXHBC_H


#include "Element.h"

class MatrixHBC {
public:
    double alpha;
    double Px[8], Py[8];
    double matrixH[4][4];

    MatrixHBC();
    ~MatrixHBC();

    void CalculateMatrixHBC(Element element);

};

#endif //FINITEELEMENTMETHOD_MATRIXHBC_H
