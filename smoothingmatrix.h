#ifndef SMOOTHINGMATRIX_H
#define SMOOTHINGMATRIX_H
#include <QVector3D>
#include <Eigen/Core>

using namespace std;

class SmoothingMatrix
{
private:
    Eigen::MatrixXf weightMatrix;
    vector<float> vertexMatrix;
    Eigen::MatrixXf resultMatrix;

public:
    SmoothingMatrix(int rows, int cols);
    void addWeightRow(int row, int startCol, vector<float> val);
    void addVertex(int row, int col, float x, float y, float z);
};

#endif // SMOOTHINGMATRIX_H
