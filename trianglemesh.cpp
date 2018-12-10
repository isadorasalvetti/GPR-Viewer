#include <iostream>
#include <math.h>
#include "trianglemesh.h"


using namespace std;

const float PI = 3.1415f;

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

float cosineToCotangent(float cos){
   return 1/tan(acos(cos));
}

//*****************************************
//Init / load and destroy
//*****************************************

TriangleMesh::TriangleMesh() : vboVertices(QOpenGLBuffer::VertexBuffer),
										 vboNormals(QOpenGLBuffer::VertexBuffer),
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
                    2, 3, 7, 7, 6, 2
                  };

	int i;

	for(i=0; i<8; i+=1)
		addVertex(0.5f * QVector3D(vertices[3*i], vertices[3*i+1], vertices[3*i+2]));
	for(i=0; i<12; i++)
		addTriangle(faces[3*i], faces[3*i+1], faces[3*i+2]);
}

bool TriangleMesh::init(QOpenGLShaderProgram *program)
{
    vector<QVector3D> replicatedVertices, normals, colors;
	vector<unsigned int> perFaceTriangles;

    colors.resize(vertices.size(), QVector3D(0.5, 0.5, 0.5));

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

	eboTriangles.destroy();
	eboTriangles.create();
	if(eboTriangles.isCreated())
		eboTriangles.bind();
	else
		return false;
	eboTriangles.setUsagePattern(QOpenGLBuffer::StaticDraw);

    fillVBOs(replicatedVertices, normals, perFaceTriangles);

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

    vector<QVector3D> colors;
    vector<float> curvatures;

    curvatures = GaussianCurvature();
    GetColors(colors, curvatures);
    colors = buildReplicatedColors(colors);
    updateColors(colors);
}

void TriangleMesh::DisplayMeanCurvature(){
    vector<QVector3D> colors;
    vector<float> curvatures;

    curvatures = MeanCurvature();
    GetColors(colors, curvatures);
    colors = buildReplicatedColors(colors);
    updateColors(colors);
}

void TriangleMesh::IteractiveSmoothing(int nSteps){
    for (int i = 0; i < nSteps; i++){
        IteractiveSmoothingStep();
    }
}

void TriangleMesh::BiIteractiveSmoothing(int nSteps){
    for (int i = 0; i < nSteps; i++){
        BiIteractiveSmoothingStep();
    }
}

void TriangleMesh::GlobalSmoothing(int percent){
    buildSmoothingMatrix(percent);
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

    vboColors.bind();
    vboColors.allocate(&normals[0], 3 * sizeof(float) * normals.size());
    vboColors.release();

	eboTriangles.bind();
	eboTriangles.allocate(&perFaceTriangles[0], sizeof(int) * perFaceTriangles.size());
	eboTriangles.release();
}

void TriangleMesh::updateColors(vector<QVector3D> &newColors)
{
    vboColors.bind();
    vboColors.allocate(&newColors[0], 3 * sizeof(float) * newColors.size());
    vboColors.release();
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

vector<int> TriangleMesh::GetVertexNeighboors(unsigned int vert){
    vector<int> neighboors;

    //get corner of this vertex
    int initialCorner = next(cornerVertex[vert]);
    int nextCorner = initialCorner;

    //find vertices surrounding this corner
    while(true){
        neighboors.push_back(nextCorner);
        nextCorner = previous(cornersTable[nextCorner]);
        if (nextCorner == -1) {break; qDebug("Non manifold mesh!");} //non manifold mesh
        if (nextCorner == initialCorner) break; //corner already in neighboors, end of loop
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
        //qDebug()<<"ggrrrargahrj";
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
        vector<int> neighboorhood = GetVertexNeighboors(v);
        for (int i = 0; i < neighboorhood.size(); i++){
            QVector3D vn = vertices[v].normalized(); QVector3D vi = vertices[triangles[neighboorhood[i]]].normalized();
            float weight = QVector3D::dotProduct(vn, vi); //TODO!
            Lv += weight*(vertices[triangles[neighboorhood[i]]]-vertices[v]);
        }
        return Lv;
    }
}

void TriangleMesh::IteractiveSmoothingStep(){
    vector<QVector3D> newVertices;
    float g = 0.3; // g E [0, 0.7]
    newVertices.resize(vertices.size());
    for (int i = 0; i < vertices.size(); i++){
        QVector3D Lv = ComputeLaplacian(i, true);
        newVertices[i] = vertices[i] + g*Lv;
    }
    vertices = newVertices;
    updateVertices();
}

void TriangleMesh::BiIteractiveSmoothingStep(){
    vector<QVector3D> newVertices;
    float g = 0.3; // g E [0, 0.7]
    newVertices.resize(vertices.size());
    for (int i = 0; i < vertices.size(); i++){
        vertices[i] = vertices[i] + g*ComputeLaplacian(i, true);
        QVector3D Lv2 = ComputeLaplacian(i, true);
        newVertices[i] = vertices[i] - g*Lv2;
    }
    vertices = newVertices;
    updateVertices();
}

//*****************************************
// L3 - GlobalSmoothing
//*****************************************

void TriangleMesh::buildSmoothingMatrix(int validHeight){
    SmoothingMatrix myMatrix(vertices.size());
    vector<bool> vertIsVariable(vertices.size(), false);
    vector<float> rowWeights;
    vector<int> rowId;

    //Select variable vertices
    for (int i = 0; i < vertices.size(); i++){
        if (vertices[i][1] > validHeight) vertIsVariable[i] = true; //variable condition
    }

    //Add weigths to matrix
    for (int i = 0; i < vertices.size(); i++){
        vector<int> neighboorhood = GetVertexNeighboors(i);
        float weight = 1/(float)neighboorhood.size();

        //Create row of weights
        rowWeights.push_back(-1);
        rowId.push_back(i);
        for (int j = 0; j < neighboorhood.size(); j++){
            rowWeights.push_back(weight);
            rowId.push_back(triangles[neighboorhood[j]]);
        }

        //add weight row and vertex to matrices
        myMatrix.addWeightRow(rowId, rowWeights);
        myMatrix.addVertex(vertices[i], i);
        rowWeights.clear();
        rowId.clear();
    }
    cout << "Smoothing matrix generated" << endl;
    vertices = myMatrix.solve();
    updateVertices();
}

//*****************************************
// L4 - Noise Magnification
//*****************************************


//*****************************************
// L5 - Discrete Harmonic Map
//*****************************************


//*****************************************
//Colors
//*****************************************

void TriangleMesh::GetColors(vector<QVector3D> &vertColors, vector<float>&vertCurvatures){
    vertColors.resize(vertCurvatures.size());
    float minCurv = *std::min_element(vertCurvatures.cbegin(), vertCurvatures.cend());
    float maxCurv = *std::max_element(vertCurvatures.cbegin(), vertCurvatures.cend());
    float maxi = std::max(abs(minCurv), abs(maxCurv));
    for (unsigned int i = 0; i<vertColors.size(); i++){
        if (vertCurvatures[i] > 0) vertColors[i] = QVector3D(sqrt(vertCurvatures[i]/maxi), 0, 0);
        else if (vertCurvatures[i] < 0) vertColors[i] = QVector3D(0, sqrt(-vertCurvatures[i]/maxi), 0);
        else vertColors[i] = QVector3D(0, 0, 0);
    }
}



