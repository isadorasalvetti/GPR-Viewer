#include "smoothingmatrix_s.h"
/*
SparseSmoothingMatrix::SparseSmoothingMatrix(int vertSize){
    weightMatrix.resize(0,1);
    resultXMatrix.resize(vertSize);
    resultYMatrix.resize(vertSize);
    resultZMatrix.resize(vertSize);
}

void SparseSmoothingMatrix::addWeightRow(vector<TEntry> entry){

}

void SparseSmoothingMatrix::buildMatrix(){

}

void SparseSmoothingMatrix::addVertex(QVector3D vertex, int i){
    resultXMatrix[i] = vertex[0];
    resultYMatrix[i] = vertex[1];
    resultZMatrix[i] = vertex[2];
}

vector<QVector3D> SparseSmoothingMatrix::solve(){
    //Lefthand side
    Eigen::MatrixXf weightMatrixT = weightMatrix.transpose(); //AtA
    Eigen::LLT<Eigen::MatrixXf> lltweightMatrix(weightMatrixT * weightMatrix);

    //Solve for X
    Eigen::VectorXf xCoords = lltweightMatrix.solve(weightMatrixT * resultXMatrix);
    //Solve for Y
    Eigen::VectorXf yCoords = lltweightMatrix.solve(weightMatrixT * resultYMatrix);
    //Solve for Z
    Eigen::VectorXf zCoords = lltweightMatrix.solve(weightMatrixT * resultZMatrix);

    vector<QVector3D> returnVector(xCoords.size());
    for (unsigned int i = 0; i < returnVector.size(); i++) returnVector[i] = QVector3D(xCoords[i], yCoords[i], zCoords[i]);

    return returnVector;
}
*/
