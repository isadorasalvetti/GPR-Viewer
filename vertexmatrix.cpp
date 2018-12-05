#include "vertexmatrix.h"

VertexMatrix::VertexMatrix(int rows, int cols){
    vertexMatrix = Eigen::MatrixXf(rows, cols);
}

VertexMatrix::addValue(int row, int col, float val){
    vertexMatrix[row][col] = val;
}
