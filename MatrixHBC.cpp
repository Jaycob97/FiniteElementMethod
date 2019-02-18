#include "stdafx.h"
#include <cmath>
#include <iostream>
#include "MatrixHBC.h"

MatrixHBC::MatrixHBC() {
        Px[0] = -1/sqrt(3);
        Py[0] = -1;

        Px[1] = 1/sqrt(3);
        Py[1] = -1;

        Px[2] = 1;
        Py[2] = -1/sqrt(3);

        Px[3] = 1;
        Py[3] = 1/sqrt(3);

        Px[4] = 1/sqrt(3);
        Py[4] = 1;

        Px[5] = -1/sqrt(3);
        Py[5] = 1;

        Px[6] = -1;
        Py[6] = 1/sqrt(3);

        Px[7] = -1;
        Py[7] = -1/sqrt(3);

        for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                       matrixH[i][j] = 0;
                }
        }
}

MatrixHBC::~MatrixHBC() {

}

void MatrixHBC::CalculateMatrixHBC(Element element) {

        double sum[4][4][4], temp1[4], temp2[4], length[4], detJ[4];

        length[0] = sqrt(pow(element.nodeID[1].x-element.nodeID[0].x, 2) + pow(element.nodeID[1].y - element.nodeID[0].y, 2));
        length[1] = sqrt(pow(element.nodeID[1].x-element.nodeID[2].x, 2) + pow(element.nodeID[1].y - element.nodeID[2].y, 2));
        length[2] = sqrt(pow(element.nodeID[2].x-element.nodeID[3].x, 2) + pow(element.nodeID[2].y - element.nodeID[3].y, 2));
        length[3] = sqrt(pow(element.nodeID[0].x-element.nodeID[3].x, 2) + pow(element.nodeID[0].y - element.nodeID[3].y, 2));

        detJ[0] = length[0]/2;
        detJ[1] = length[1]/2;
        detJ[2] = length[2]/2;
        detJ[3] = length[3]/2;

        for(int i=0; i<4; i++){
                temp1[0] = (0.25 * (1 - Px[2*i]) * (1 - Py[2*i]));
                temp1[1] = (0.25 * (1 + Px[2*i]) * (1 - Py[2*i]));
                temp1[2] = (0.25 * (1 + Px[2*i]) * (1 + Py[2*i]));
                temp1[3] = (0.25 * (1 - Px[2*i]) * (1 + Py[2*i]));

                temp2[0] = (0.25 * (1 - Px[2*i+1]) * (1 - Py[2*i+1]));
                temp2[1] = (0.25 * (1 + Px[2*i+1]) * (1 - Py[2*i+1]));
                temp2[2] = (0.25 * (1 + Px[2*i+1]) * (1 + Py[2*i+1]));
                temp2[3] = (0.25 * (1 - Px[2*i+1]) * (1 + Py[2*i+1]));

                for(int j=0;j<4;j++){
                        for(int k=0;k<4;k++){
                                sum[i][j][k] = ((alpha * temp1[j] * temp1[k]) + (alpha * temp2[j] * temp2[k])) * detJ[i];
                        }
                }
        }
        for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                        for(int k=0;k<4;k++){
                                matrixH[i][j] += element.isSurfaceHeated[k] * sum[k][i][j];
                        }
                }
        }
}

