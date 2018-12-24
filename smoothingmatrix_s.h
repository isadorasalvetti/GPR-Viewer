#ifndef SMOOTHINGMATRIX_S_H
#define SMOOTHINGMATRIX_S_H
#include <QVector3D>
#include <QVector>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

using namespace std;
typedef Eigen::Triplet<double> TEntry;

class SparseSmoothingMatrix
{/*
private:
    vector<TEntry> weightEntries;
    Eigen::SparseMatrix<float> weightMatrix;
    vector<float> vertexMatrix;
    Eigen::VectorXf resultXMatrix;
    Eigen::VectorXf resultYMatrix;
    Eigen::VectorXf resultZMatrix;


public:
    SmoothingMatrix(int vertSize);
    vector<QVector3D> solve();
    void addWeightRow(vector<TEntry> entry);
    void addVertex(QVector3D vertex, int i);
    void buildMatrix();
    void addResultEntry();*/
};

#endif // SMOOTHINGMATRIX_S_H
