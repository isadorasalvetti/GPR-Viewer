#include "smoothingmatrix.h"

SmoothingMatrix::SmoothingMatrix(){
}

void SmoothingMatrix::addWeightRow(vector<int> rowId, vector<float> rowWeights){
    int cols = weightMatrix.cols();
    int rows = weightMatrix.rows();

    //Resize matrix to accomodate a new row
    weightMatrix.conservativeResize(rows+1, cols);

    //Add coeficients
    for (int i=0; i<rowId.size(); i++){
        //rowTd[0]/ rowWeights[0] is the id of the current vertex. The remaining items are its neighboors.
        weightMatrix(rowId[0], rowId[i]) = rowWeights[i];
    }
}
