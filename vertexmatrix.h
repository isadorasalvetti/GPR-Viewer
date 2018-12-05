#ifndef VERTEXMATRIX_H
#define VERTEXMATRIX_H

#include <QVector3D>
#include "Eigen/Core/Matrix.h"

using namespace std;

class VertexMatrix
{
private:
    Eigen::MatrixXf weightMatrix;
    vector<float> vertexMatrix;
    Eigen::MatrixXf resultMatrix;

public:
    VertexMatrix(int rows, int cols);
    void addWeightRow(int row, int startCol, float[5] val);
    void addVertex(int row, int col, float x, float y, float z);
};

#endif // VERTEXMATRIX_H
