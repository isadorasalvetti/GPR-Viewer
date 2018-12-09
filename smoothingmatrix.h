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
    Eigen::VectorXf resultXMatrix;
    Eigen::VectorXf resultYMatrix;
    Eigen::VectorXf resultZMatrix;


public:
    SmoothingMatrix(int vertSize);
    void addWeightRow(vector<int> rowId, vector<float> rowWeights);
    void addVertex(QVector3D vertex, int i);
    void addResultEntry();
};

#endif // SMOOTHINGMATRIX_H
