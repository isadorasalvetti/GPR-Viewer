#include <iostream>
#include <math.h>
#include "trianglemesh.h"
#include <Eigen/Sparse>


using namespace std;

const double PI = 3.14159265359;

//*****************************************
// CornerEdges
//*****************************************

int next(int corner)
{
	return 3 * (corner / 3) + (corner + 1) % 3;
}

int previous(int corner)
{
	return 3 * (corner / 3) + (corner + 2) % 3;
}

bool CornerEdge::Compare(CornerEdge e1, CornerEdge e2){
    bool result = (e1.edgeMin == e2.edgeMin && e1.edgeMax == e2.edgeMax);
    return result;
}

//*****************************************
//Utilities
//*****************************************

pair<QVector3D, QVector3D> getBoundingBox(vector<QVector3D> &verts){
    float vertCount = verts.size();
    float x0 = 0; float y0 = 0; float z0 = 0; float x1 = 0; float y1 = 0; float z1 = 0;
    for (int i =0; i<vertCount; i++){
        x0 = min(x0, verts[i][0]);
        y0 = min(y0, verts[i][1]);
        z0 = min(z0, verts[i][2]);
        x1 = max(x1, verts[i][0]);
        y1 = max(y1, verts[i][1]);
        z1 = max(z1, verts[i][2]);
    }
    return pair<QVector3D, QVector3D>(QVector3D(x0, y0, z0), QVector3D(x1, y1, z1));
}

pair<float, float> getBoundingBoxHeight(vector<QVector3D> &verts){
    float vertCount = verts.size();
    float x0 = 0; float y0 = 0; float z0 = 0; float x1 = 0; float y1 = 0; float z1 = 0;
    for (int i =0; i<vertCount; i++){
        y0 = min(y0, verts[i][1]);
        y1 = max(y1, verts[i][1]);
    }
    return pair<float, float>(y0, y1);
}

const float eps = 0.000001;
const float cotan_max = cos( eps ) / sin( eps );

float cosineToCotangent(float cos){
   float bt = sqrt(1.0f - cos * cos);
   //if (bt < 0.0001f and bt > -0.0001f) return 0;
   return abs(cos / bt);
}

float getCotangent(QVector3D va, QVector3D vb){
    va.normalize(); vb.normalize();
    float sin = (QVector3D::crossProduct(va, vb)).length();
    if (abs(sin) < 0.000001) return 0;
    float cos = QVector3D::dotProduct(va, vb);

    float cot = cos/sin;
    if (cot < 0) cot = -cot;
    cot = max(cot, -cotan_max);
    return cot;
}

//*****************************************
//Init / load and destroy
//*****************************************

TriangleMesh::TriangleMesh() : vboVertices(QOpenGLBuffer::VertexBuffer),
                               vboNormals(QOpenGLBuffer::VertexBuffer),
                               vboColors(QOpenGLBuffer::VertexBuffer),
                               vboUV(QOpenGLBuffer::VertexBuffer),
                               eboTriangles(QOpenGLBuffer::IndexBuffer)
{
}


void TriangleMesh::addVertex(const QVector3D &position)
{
	vertices.push_back(position);
}

void TriangleMesh::addTriangle(int v0, int v1, int v2)
{
	triangles.push_back(v0);
	triangles.push_back(v1);
	triangles.push_back(v2);
}

void TriangleMesh::buildCube()
{
	GLfloat vertices[] = {-1, -1, -1,
                          1, -1, -1,
                          1,  1, -1,
                         -1,  1, -1,
                         -1, -1,  1,
                          1, -1,  1,
                          1,  1,  1,
                         -1,  1,  1
                        };

	GLuint faces[] = {3, 1, 0, 3, 2, 1,
                    5, 6, 7, 4, 5, 7,
                    7, 3, 0, 0, 4, 7,
                    1, 2, 6, 6, 5, 1,
                    0, 1, 4, 5, 4, 1,
//                    2, 3, 7, 7, 6, 2
                  };

	int i;

	for(i=0; i<8; i+=1)
		addVertex(0.5f * QVector3D(vertices[3*i], vertices[3*i+1], vertices[3*i+2]));
    for(i=0; i<10; i++)
		addTriangle(faces[3*i], faces[3*i+1], faces[3*i+2]);
}

bool TriangleMesh::init(QOpenGLShaderProgram *program){
    pg = program;
    vector<QVector3D> replicatedVertices, normals, colors;
	vector<unsigned int> perFaceTriangles;

    colors.resize(vertices.size(), QVector3D(1, 1, 1));

    buildReplicatedVertices(replicatedVertices, normals, perFaceTriangles);
    buildCornerTable();

	program->bind();

	vao.destroy();
	vao.create();
	if(vao.isCreated())
		vao.bind();
	else
		return false;

    vboVertices.destroy();
    vboVertices.create();
    if(vboVertices.isCreated())
        vboVertices.bind();
    else
        return false;
    vboVertices.setUsagePattern(QOpenGLBuffer::StaticDraw);
    program->enableAttributeArray(0);
    program->setAttributeBuffer(0, GL_FLOAT, 0, 3, 0);

    vboNormals.destroy();
    vboNormals.create();
    if(vboNormals.isCreated())
        vboNormals.bind();
    else
        return false;
    vboNormals.setUsagePattern(QOpenGLBuffer::StaticDraw);
    program->enableAttributeArray(1);
    program->setAttributeBuffer(1, GL_FLOAT, 0, 3, 0);

    vboColors.destroy();
    vboColors.create();
    if(vboColors.isCreated())
        vboColors.bind();
    else
        return false;
    vboColors.setUsagePattern(QOpenGLBuffer::StaticDraw);
    program->enableAttributeArray(2);
    program->setAttributeBuffer(2, GL_FLOAT, 0, 3, 0);

    vboUV.destroy();
    vboUV.create();
    if(vboUV.isCreated())
        vboUV.bind();
    else
        return false;
    vboUV.setUsagePattern(QOpenGLBuffer::StaticDraw);
    program->enableAttributeArray(3);
    program->setAttributeBuffer(3, GL_FLOAT, 0, 2, 0);

	eboTriangles.destroy();
	eboTriangles.create();
	if(eboTriangles.isCreated())
		eboTriangles.bind();
	else
		return false;
	eboTriangles.setUsagePattern(QOpenGLBuffer::StaticDraw);

    fillVBOs(replicatedVertices, normals, perFaceTriangles);
    vector<QVector3D> repColors = buildReplicatedColors(colors);
    updateColors(repColors);
    backupVertices = vertices;

	vao.release();
	program->release();

	return true;
}

void TriangleMesh::destroy()
{
	vao.destroy();
	vboVertices.destroy();
    vboColors.destroy();
	vboNormals.destroy();
	eboTriangles.destroy();

	vertices.clear();
	triangles.clear();
}

//*****************************************
// Buttons
//*****************************************

void TriangleMesh::DisplayGaussianCurvature(){
    pg->bind();
    pg->setUniformValue("useUV", false);

    vector<QVector3D> colors;
    vector<float> curvatures;

    curvatures = GaussianCurvature();
    GetColors(colors, curvatures);
    colors = buildReplicatedColors(colors);
    updateColors(colors);
    pg->release();
}

void TriangleMesh::DisplayMeanCurvature(){
    pg->bind();
    pg->setUniformValue("useUV", false);

    vector<QVector3D> colors;
    vector<float> curvatures;

    curvatures = MeanCurvature();
    GetColors(colors, curvatures);
    colors = buildReplicatedColors(colors);
    updateColors(colors);
    pg->release();
}

void TriangleMesh::IteractiveSmoothing(int nSteps, bool uw){
    for (int i = 0; i < nSteps; i++){
        IteractiveSmoothingStep(uw);
    }
}

void TriangleMesh::BiIteractiveSmoothing(int nSteps, bool uw){
    for (int i = 0; i < nSteps; i++){
        BiIteractiveSmoothingStep(uw);
    }
}

void TriangleMesh::GlobalSmoothing(float percent, bool type){
    pair<float, float> bBoxHt = getBoundingBoxHeight(vertices);
    float diff = bBoxHt.second - bBoxHt.first;
    float cutOffHeight = bBoxHt.first + diff*percent;
    if (type) solveSparseSmoothing(getCutOffVerticesPrcnt(percent));
    else solveSparseSmoothing(getCutOffVertices(cutOffHeight));
    updateVertices();
}

void TriangleMesh::DetMagnification(QVector3D l){
    getNoise(l[0], l[1], l[2]);
}
//**********
vector<bool> convertToBool(vector<int>constVerts, int size){
    vector<bool> vertIsVariable(size, true);
    for (int i =0; i<constVerts.size(); i++)
        vertIsVariable[constVerts[i]] = false;
    return vertIsVariable;
}
//*************

void TriangleMesh::Parametrization(){
    pg->bind();
    pg->setUniformValue("useUV", true);
    if (!uvsComputed){
        vector<QVector3D> saveVertices = vertices;
        vector<int> borderVertices = getBoundary();
        CreateMapBorder(borderVertices);
        solveSparseSmoothing(convertToBool(borderVertices, vertices.size()));
        uvsComputed = true;
        uvs = vertices;
        vertices = backupVertices;
    }
    vector<QVector2D> newUVs = buildReplicatedUVs(uvs);
    updateUVs(newUVs);
    pg->release();
}

void TriangleMesh::DisplayParametrization(){
    pg->bind();
    pg->setUniformValue("useUV", true);
    if (!uvsComputed){
        vector<int> borderVertices = getBoundary();
        CreateMapBorder(borderVertices);
        solveSparseSmoothing(convertToBool(borderVertices, vertices.size()));
        uvsComputed = true;
        uvs = vertices;
    }
    else vertices = uvs;
    updateVertices();
    vector<QVector2D> newUVs = buildReplicatedUVs(vertices);
    updateUVs(newUVs);
    pg->release();
}


void TriangleMesh::Reset(){
    vertices = backupVertices;
    uvsComputed = false;
    updateVertices();
}

//*****************************************
// Render
//*****************************************

void TriangleMesh::render(QOpenGLFunctions &gl)
{
	vao.bind();
	eboTriangles.bind();
    gl.glDrawElements(GL_TRIANGLES, triangles.size(), GL_UNSIGNED_INT, nullptr);
	vao.release();
}

//****************************************
// Building and filling buffers
//*****************************************

vector<QVector3D> TriangleMesh::buildReplicatedColors(vector<QVector3D> currColors){
    vector<QVector3D> repColors;
    for(unsigned int i=0; i<triangles.size(); i+=3)    {
        repColors.push_back(currColors[triangles[i]]);
        repColors.push_back(currColors[triangles[i+1]]);
        repColors.push_back(currColors[triangles[i+2]]);
    }
    return repColors;
}

vector<QVector2D> TriangleMesh::buildReplicatedUVs(vector<QVector3D> currUVs){
    vector<QVector2D> repUVs;
    for(unsigned int i=0; i<triangles.size(); i+=3)    {
        repUVs.push_back(
        QVector2D(currUVs[triangles[i]][0], currUVs[triangles[i]][2]));

        repUVs.push_back(
        QVector2D(currUVs[triangles[i+1]][0], currUVs[triangles[i+1]][2]));

        repUVs.push_back(
        QVector2D(currUVs[triangles[i+2]][0], currUVs[triangles[i+2]][2]));
    }
    return repUVs;
}

void TriangleMesh::buildReplicatedVertices(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles)
{
	normals.resize(triangles.size());

	for(unsigned int i=0; i<triangles.size(); i+=3)
	{
		replicatedVertices.push_back(vertices[triangles[i]]);
		replicatedVertices.push_back(vertices[triangles[i+1]]);
		replicatedVertices.push_back(vertices[triangles[i+2]]);

		QVector3D N = QVector3D::crossProduct(vertices[triangles[i+1]] - vertices[triangles[i]], vertices[triangles[i+2]] - vertices[triangles[i]]);
		N.normalize();
		normals[i] = N;
		normals[i+1] = N;
		normals[i+2] = N;

		perFaceTriangles.push_back(perFaceTriangles.size());
		perFaceTriangles.push_back(perFaceTriangles.size());
		perFaceTriangles.push_back(perFaceTriangles.size());
	}
}

void TriangleMesh::fillVBOs(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles)
{
	vboVertices.bind();
	vboVertices.allocate(&replicatedVertices[0], 3 * sizeof(float) * replicatedVertices.size());
	vboVertices.release();

    vboNormals.bind();
    vboNormals.allocate(&normals[0], 3 * sizeof(float) * normals.size());
    vboNormals.release();

	eboTriangles.bind();
	eboTriangles.allocate(&perFaceTriangles[0], sizeof(int) * perFaceTriangles.size());
	eboTriangles.release();
}

void TriangleMesh::updateColors(vector<QVector3D> &newColors){
    vboColors.bind();
    vboColors.allocate(&newColors[0], 3 * sizeof(float) * newColors.size());
    vboColors.release();
}

void TriangleMesh::updateUVs(vector<QVector2D> &newUVs){
    vboUV.bind();
    vboUV.allocate(&newUVs[0], 2 * sizeof(float) * newUVs.size());
    vboUV.release();
}

void TriangleMesh::updateVertices() {
    vector<QVector3D> replicatedVertices, normals;
    vector<unsigned int> perFaceTriangles;
    buildReplicatedVertices(replicatedVertices, normals, perFaceTriangles);
    fillVBOs(replicatedVertices, normals, perFaceTriangles);
}


//*****************************************
//Corner table and neighborhood
//*****************************************

void TriangleMesh::buildCornerTable(){

    // Get/ store edge-corner, vertex-corner data.
    unsigned int NoTris = static_cast<unsigned int>(triangles.size()/3);
    cornerVertex.resize(vertices.size());

    vector<CornerEdge> allCornersEdges;
    for (unsigned int i = 0; i < NoTris; i++){ //(1)
        //Triangle vertices
        int v0 = triangles[i*3 + 0];
        int v1 = triangles[i*3 + 1];
        int v2 = triangles[i*3 + 2];

        //Edges
        //i
        int e1Min = std::min(v0, v1);
        int e1Max = std::max(v0, v1);
        //i + 1
        int e2Min = std::min(v1, v2);
        int e2Max = std::max(v1, v2);
        //i + 2
        int e3Min = std::min(v0, v2);
        int e3Max = std::max(v0, v2);

        //Assign Corner to Edges
        allCornersEdges.push_back(CornerEdge(e2Min, e2Max, i*3 + 0));
        allCornersEdges.push_back(CornerEdge(e3Min, e3Max, i*3 + 1));
        allCornersEdges.push_back(CornerEdge(e1Min, e1Max, i*3 + 2));

        //Assign Coner from Vertex
        cornerVertex[v0]=i*3+0;
        cornerVertex[v1]=i*3+1;
        cornerVertex[v2]=i*3+2;
        }

    //Build Corner table
    cornersTable.resize(allCornersEdges.size());
    for (int i =0; i < allCornersEdges.size(); i++){
        CornerEdge e1 = allCornersEdges[i];
        cornersTable[i] = -1;

        for (int j =0; j < allCornersEdges.size(); j++){
            CornerEdge e2 = allCornersEdges[j];

            if (i!=j && CornerEdge::Compare(e1, e2)){
                cornersTable[i] = e2.cornerVertex;
                cornersTable[j] = e1.cornerVertex;
                break;
                }
            }
        }
    std::cout<<"Corner table generated"<<std::endl;
}

vector<int> TriangleMesh::getBoundary(){
    list<int> singleCorners;
    isInBorder.resize(vertices.size(), false);
    for (int i=0; i< cornersTable.size(); i++){
        if (cornersTable[i] == -1){ //found first border vertex.
            singleCorners.push_back(i);
            isInBorder[i] = true;
        }
    }

    if (singleCorners.size()==0) {
        cout<<"Error in parametrization: Mesh must be not watertight";
        exception e; throw e;
    }

    vector<int> borderEdge;

    int faceI = singleCorners.front();
    singleCorners.erase(singleCorners.begin());
    int vertA = triangles[next(faceI)];
    int vertB = triangles[next(next(faceI))];
    borderEdge.push_back(vertA);
    borderEdge.push_back(vertB);

    while(singleCorners.size() > 0)
    for (auto j = singleCorners.begin();
         j != singleCorners.end();
         j++){
            int faceJ = *j;
            vertA = triangles[next(faceJ)];
            if (vertA == vertB) {
                vertB = triangles[next(next(faceJ))];
                borderEdge.push_back(vertB);
                singleCorners.erase(j);
                break;
            }
        }

//    //Color verification:
//    vector<QVector3D> colors (vertices.size(), QVector3D(0, 0, 0));
//    for (int i =0; i< borderEdge.size(); i++){
//        colors[borderEdge[i]] = QVector3D(1, 0, 0);
//        vector<QVector3D> rColors = buildReplicatedColors(colors);
//        updateColors(rColors);
//        gl->update();
//    }
    return borderEdge;
}

vector<int> TriangleMesh::GetVertexNeighboors(unsigned int vert){
    vector<int> neighboors;
    vector<int> temp;

    //get corner of this vertex
    int nextCorner = next(cornerVertex[vert]);
    neighboors.push_back(nextCorner);

    //find vertices surrounding this corner
    while(true){
        nextCorner = next(nextCorner);
        nextCorner = cornersTable[nextCorner];
        if (nextCorner < 0) { //Mesh is non-manifold. Must check border
            nextCorner = previous(cornerVertex[vert]);
            neighboors.push_back(nextCorner);
            while(true){
                nextCorner = previous(nextCorner);
                nextCorner = cornersTable[nextCorner];
                if (nextCorner < 0) break;
                temp.push_back(nextCorner);
            }
            for (int i = temp.size()-1; i >= 0; i--)
                neighboors.push_back(temp[i]);
            break;
        }
    if (nextCorner == next(cornerVertex[vert])) break; //corner already in neighboors, end of loop
    neighboors.push_back(nextCorner);
    }
    return neighboors;
}

//*****************************************
// L1 - Curvatures
//*****************************************

vector<float> TriangleMesh::GaussianCurvature(){
    vector<float> perVertCurvature(vertices.size());

    //kg = 2pi - sum(angle) / area.
    for (unsigned int i = 0; i < vertices.size(); i++){
        //get neighboorhood of the vertex
        vector<int> neighboors = GetVertexNeighboors(i);

        //get sum of neighbooring angles
        float angle = 0;

        //get area of triangle
        float area = 0;
        for (unsigned int j=0; j < neighboors.size(); j++){
            QVector3D v1 = - vertices[i] + vertices[triangles[neighboors[j]]];
            QVector3D v2 = - vertices[i] + vertices[triangles[neighboors[(j+1) % neighboors.size()]]];
            float ar = QVector3D::crossProduct(v1, v2).length()/2;
            v1.normalize(); v2.normalize();
            float a = acos(QVector3D::dotProduct(v1, v2));
            //qDebug()<<a;
            angle += a;
            area += ar;
        }
        float curvature = (2*PI - angle);
        perVertCurvature[i] = curvature;
        //qDebug()<<curvature;
    }
    return perVertCurvature;
}

vector<float> TriangleMesh::MeanCurvature(){
    //kh = len(1/2A * sum(cot(alpha)+cot(beta))*(vi-vj)).
    vector<float> perVertCurvature(vertices.size());

    for (unsigned int i = 0; i < vertices.size(); i++){
        vector<int> neighboors = GetVertexNeighboors(i);
        QVector3D sum = QVector3D(0,0,0);
        float area = 0;
        QVector3D v0 = vertices[i];
        for (unsigned int j=0; j < neighboors.size(); j++){
            QVector3D vmin = vertices[triangles[neighboors[j]]];
            QVector3D vi = vertices[triangles[neighboors[(j+1) % neighboors.size()]]];
            QVector3D vplus = vertices[triangles[neighboors[(j+2) % neighboors.size()]]];

            QVector3D iEdge = vplus - v0; QVector3D iEdge2 = vplus - vi;
            QVector3D jEdge = vmin - v0; QVector3D jEdge2 = vmin - vi;

            iEdge.normalize(); jEdge.normalize();
            iEdge2.normalize(); jEdge2.normalize();

            float ar = QVector3D::crossProduct(v0, vmin).length()/2;

            float alpha = QVector3D::dotProduct(iEdge, iEdge2);
            float beta = QVector3D::dotProduct(jEdge, jEdge2);
            sum += (cosineToCotangent(alpha)+cosineToCotangent(beta)) * (v0-vi);
            area += ar;
        }
        perVertCurvature[i] = (1/(2*area)) * sum.length();
    }
    return perVertCurvature;

}

//*****************************************
// L2 - Smoothing
//*****************************************

int getNext(const int &i, const int &count){
    int nxt = i+1;
    if (nxt < count) return nxt;
    else return 0;
}

int getPrev(const int &i, const int &count){
    int prev = i-1;
    if (prev >= 0) return prev;
    else return count-1;
}

QVector3D TriangleMesh::ComputeLaplacian(int v, bool uniform){
    //Returns Laplacian of vertex v
    QVector3D Lv = QVector3D(0, 0, 0);
    if (uniform){
        vector<int> neighboorhood = GetVertexNeighboors(v);
        float weight = 1/(float)neighboorhood.size();
        for (int i = 0; i < neighboorhood.size(); i++) Lv += weight*(vertices[triangles[neighboorhood[i]]]-vertices[v]);
        return Lv;
    }

    else { //CotangentWeights
        vector<int> neighborhood = GetVertexNeighboors(v);
        int nSize = neighborhood.size();
        for (int i = 0; i < nSize; i++){
            Lv += getCotangentWeight(v, i, neighborhood)*(vertices[triangles[neighborhood[i]]]-vertices[v]);
        }
        return Lv;
    }
}

float TriangleMesh::getCotangentWeight(int v, int i, vector<int> &neighborhood){
    int nSize = neighborhood.size();
    QVector3D vn = vertices[v]; QVector3D vi = vertices[triangles[neighborhood[i]]];
    int indexKa = triangles[neighborhood[getPrev(i, nSize)]];
    int indexKb = triangles[neighborhood[getNext(i, nSize)]];

    QVector3D ka = vertices[indexKa];
    QVector3D kb = vertices[indexKb];
    QVector3D vnka = (ka - vn); QVector3D vnkb = (kb - vn);
    QVector3D vika = (ka - vi); QVector3D vikb = (kb - vi);
    float w1 = getCotangent(vnka, vika);
    float w2 = getCotangent(vnkb, vikb);

    //cout << "Weights " << w1 << ", " << w2 << endl;

   return (w1 + w2)/2;
}

void TriangleMesh::IteractiveSmoothingStep(bool weight){
    vector<QVector3D> newVertices;
    float g = 0.3; // g E [0, 0.7]
    newVertices.resize(vertices.size());
    for (int i = 0; i < vertices.size(); i++){
        QVector3D Lv = ComputeLaplacian(i, weight);
        newVertices[i] = vertices[i] + g*Lv;
    }
    vertices = newVertices;
    updateVertices();
}

void TriangleMesh::BiIteractiveSmoothingStep(bool weight){
    vector<QVector3D> newVertices;
    float g = 0.3; // g E [0, 0.7]
    newVertices.resize(vertices.size());
    for (int i = 0; i < vertices.size(); i++){
        vertices[i] = vertices[i] + g*ComputeLaplacian(i, weight);
        QVector3D Lv2 = ComputeLaplacian(i, weight);
        newVertices[i] = vertices[i] - g*Lv2;
    }
    vertices = newVertices;
    updateVertices();
}

//*****************************************
// L3 - GlobalSmoothing
//*****************************************


//vector<bool> TriangleMesh::getCutOffVertices(float validHeight){
//    //Select variable vertices
//    vector<bool> vertIsVariable(vertices.size(), false);
//    for (int i = 0; i < vertices.size(); i++){
//        vertIsVariable[i] = true or vertices[i][1] > validHeight;
//    }
//    return vertIsVariable;
//}

//void TriangleMesh::buildSmoothingMatrix(QVector3D vertIsVariable){
//    SmoothingMatrix myMatrix(vertices.size());
//    vector<float> rowWeights;
//    vector<int> rowId;

//    //Add weigths to matrix
//    for (int i = 0; i < vertices.size(); i++){
//        vector<int> neighboorhood = GetVertexNeighboors(i);
//        float weight = 1.0/neighboorhood.size();

//        //Create row of weights
//        rowWeights.push_back(-1);
//        rowId.push_back(i);
//        for (int j = 0; j < neighboorhood.size(); j++){
//            rowWeights.push_back(weight);
//            rowId.push_back(triangles[neighboorhood[j]]);
//        }

//        //add weight row and vertex to matrices
//        myMatrix.addWeightRow(rowId, rowWeights);
//        myMatrix.addVertex(vertices[i], i);
//        rowWeights.clear();
//        rowId.clear();
//    }
//    cout << "Smoothing matrix generated" << endl;
//    vertices = myMatrix.solve();
//    updateVertices();
//}

vector<bool> TriangleMesh::getCutOffVertices(float validHeight){
    // build entries
    unsigned int n = vertices.size();
    vector<bool> vertIsVariable(n, false);
    for (int i = 0; i < n; i++) {
        vertIsVariable[i] = vertices[i].y() > validHeight;
    }
    cout << "Cut off height set too " << validHeight << "." << endl;
    return vertIsVariable;
}

vector<bool> TriangleMesh::getCutOffVerticesPrcnt(float percentage){
    // build entries
    unsigned int n = vertices.size();
    vector<bool> vertIsVariable(n, false);
    for (int i = 0; i < n; i++) {
        float random = rand()/(float)RAND_MAX;
        vertIsVariable[i] = random > percentage;
    }
    cout << percentage << "of vertices set to fixed." << endl;
    return vertIsVariable;
}

void TriangleMesh::solveSparseSmoothing(vector<bool> vertIsVariable) {
    vector<Eigen::Triplet<float>> coeffs;
    Eigen::VectorXf b[3];
    unsigned int n = vertices.size();

    // init b to the max size
    for (int k = 0; k < 3; ++k) b[k] = Eigen::VectorXf(n);

    // Matrix A is m1 x m2,  m2x1 unknowns,  m1x1 b
    unsigned int m1 = 0;
    unsigned int m2 = 0;
    vector<unsigned int> j2J(n, -1);
    for (int i = 0; i < n; ++i) {
        if (vertIsVariable[i]) {
            j2J[i] = m2;
            ++m2;
        }
    }

    unsigned int I = 0;
    for (unsigned int i = 0; i < n; ++i) {
        QVector3D &v_i  = vertices[i];
        vector<int> neighbors = GetVertexNeighboors(i);
        float sumWeight_i = 0;
        float b_i[3] = {0,0,0};
        bool isRowNull = true;
        for (int l = 0; l < neighbors.size(); l++) {
            int j = triangles[neighbors[l]];
            QVector3D &v_j = vertices[j];

            float weight_i_j = 1;

            //float weight_i_j = getCotangentWeight(i, l, neighbors);
            sumWeight_i += weight_i_j;
            if (vertIsVariable[j]) {
                isRowNull = false;
                coeffs.push_back(Eigen::Triplet<float> (I, j2J[j], weight_i_j));
            } else {
                for (int k = 0; k < 3; ++k)
                    b_i[k] += -weight_i_j * v_j[k];
            }
        }
        if (vertIsVariable[i]) {
            isRowNull = false;
            coeffs.push_back(Eigen::Triplet<float> (I, j2J[i], -sumWeight_i));
        } else {
            for (int k = 0; k < 3; ++k)
                b_i[k] += sumWeight_i * v_i[k];
        }

        if (not isRowNull) {
            for (int k = 0; k < 3; ++k) b[k][I] = b_i[k];
            ++I;
        }
    }
    m1 = I;
    // solve system
    Eigen::SparseMatrix<float> A(m1,m2);
    A.setFromTriplets(coeffs.begin(), coeffs.end());

    Eigen::SparseMatrix<float> At = A.transpose();
    // Solving:
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> chol;
    chol.compute(At*A);

    if (chol.info() != Eigen::Success) {
        cout << "SOLVER PROBLEM " << chol.info() << '\n'
             << Eigen::Success << ' ' << Eigen::NumericalIssue << ' '
             << Eigen::NoConvergence << ' ' << Eigen::InvalidInput << endl;
    }
    cout << "CHOLESKY OK" << endl;


    Eigen::VectorXf sol[3];
    for (int k = 0; k < 3; ++k)
        sol[k] = chol.solve(At*b[k].block(0,0,m1,1));

    // change representation
    for (unsigned int i = 0; i < vertices.size(); ++i) {
        if (vertIsVariable[i])
            vertices[i] = QVector3D(sol[0][j2J[i]], sol[1][j2J[i]], sol[2][j2J[i]]);
    }
}

//*****************************************
// L4 - Detail Magnification
//*****************************************

void TriangleMesh::getNoise(float l1, float l2, float l3){
    int aSmthng = 10;
    vector<QVector3D> smoothedVertices = SmoothingSteps(vertices, aSmthng);

    //M = M + l1*(S0-S1) + l2*(S1-S2) + l3*(S2-S3)
    if (l1 > 0)
        for (int i = 0; i < vertices.size(); i++) vertices[i] = vertices[i] + l1*(vertices[i] - smoothedVertices[i]);
    if (l2 > 0 || l3 > 0){
        vector<QVector3D> smoothedVertices1 = SmoothingSteps(smoothedVertices, aSmthng);
        vector<QVector3D> smoothedVertices2 = SmoothingSteps(smoothedVertices1, aSmthng);
        if (l2 > 0) for (int i = 0; i < vertices.size(); i++) vertices[i] = vertices[i] + l2*(smoothedVertices[i] - smoothedVertices1[i]);
        if (l3 > 0) for (int i = 0; i < vertices.size(); i++) vertices[i] = vertices[i] + l3*(smoothedVertices1[i] - smoothedVertices2[i]);
    }
    updateVertices();
}

vector<QVector3D> TriangleMesh::SmoothingSteps(vector<QVector3D> &toSmooth, int n){
    vector<QVector3D> newVertices = toSmooth;
    for (int i = 0; i < n; i++){
        float g = 0.3; // g E [0, 0.7]
        newVertices.resize(toSmooth.size());
        for (int i = 0; i < toSmooth.size(); i++){
            QVector3D Lv = ComputeLaplacian(i, true);
            newVertices[i] = newVertices[i] + g*Lv;
        }
    }
    return newVertices;
}

//*****************************************
// L5 - Discrete Harmonic Map
//*****************************************

void TriangleMesh::CreateMapBorder(vector<int> borderVertices){
    int size = borderVertices.size();
    double angle = (2*PI)/size;

    for (int i = 0; i < size; i++){
         QVector3D p1 = QVector3D(sin(angle*i)/2+0.5f, 0, cos(angle*i)/2+0.5f);
         vertices[borderVertices[i]] = p1;
    }
}

//*****************************************
// Colors
//*****************************************

void TriangleMesh::GetColors(vector<QVector3D> &vertColors, vector<float>&vertCurvatures){
    vertColors.resize(vertCurvatures.size());
    float minCurv = *std::min_element(vertCurvatures.cbegin(), vertCurvatures.cend());
    float maxCurv = *std::max_element(vertCurvatures.cbegin(), vertCurvatures.cend());
    float maxi = std::max(abs(minCurv), abs(maxCurv));
    for (unsigned int i = 0; i<vertColors.size(); i++){
        if (vertCurvatures[i] > 0) vertColors[i] = QVector3D(sqrt(sqrt(vertCurvatures[i]/maxi)), 0, 0);
        else if (vertCurvatures[i] < 0) vertColors[i] = QVector3D(0, sqrt(sqrt(-vertCurvatures[i]/maxi)), 0);
        else vertColors[i] = QVector3D(0, 0, 0);
    }
}



