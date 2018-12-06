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
    SmoothingMatrix();
    void addWeightRow(vector<int> rowId, vector<float> rowWeights);
    void addVertex(int row, int col, float x, float y, float z);
    void addResultEntry();
};

#endif // SMOOTHINGMATRIX_H
