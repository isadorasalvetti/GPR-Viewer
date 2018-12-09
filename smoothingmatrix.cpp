#include "smoothingmatrix.h"

SmoothingMatrix::SmoothingMatrix(int vertSize){
    weightMatrix.resize(0,1);
    resultXMatrix.resize(vertSize);
    resultYMatrix.resize(vertSize);
    resultZMatrix.resize(vertSize);
}

void SmoothingMatrix::addWeightRow(vector<int> rowId, vector<float> rowWeights){
    int cols = weightMatrix.cols();
    int rows = weightMatrix.rows();

    //Resize matrix to accomodate a new row
    weightMatrix.conservativeResizeLike(Eigen::MatrixXf::Zero(rows+1, cols));
    rows += 1;

    //Add coeficients
    for (int i=0; i<rowId.size(); i++){
        //rowTd[0]/ rowWeights[0] is the id of the current vertex. The remaining items are its neighboors.
        if (rowId[i] >= weightMatrix.cols()) weightMatrix.conservativeResizeLike(Eigen::MatrixXf::Zero(rows, rowId[i]+1));
        weightMatrix(rowId[0], rowId[i]) = rowWeights[i];
    }
    int stop = 0;
}

void SmoothingMatrix::addVertex(QVector3D vertex, int i){
    resultXMatrix[i] = vertex[0];
    resultYMatrix[i] = vertex[1];
    resultZMatrix[i] = vertex[2];
}

vector<QVector3D> SmoothingMatrix::solve(){
    //Solve for X
    //Solve for Y
    //Solve for Z
}
