 #include <iostream>
#include <math.h>
#include "trianglemesh.h"


using namespace std;

const float PI = 3.1415;

int next(int corner)
{
	return 3 * (corner / 3) + (corner + 1) % 3;
}

int previous(int corner)
{
	return 3 * (corner / 3) + (corner + 2) % 3;
}

bool CornerEdge::Compare(CornerEdge e1, CornerEdge e2){
    bool result = e1.edgeMin == e2.edgeMin && e1.edgeMax == e2.edgeMax;
    return result;
}

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

    buildReplicatedVertices(replicatedVertices, normals, perFaceTriangles);

    buildCornerTable();
    float min, max;
    GetColors(colors, GaussianCurvature(min, max), min, max);

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

    vboColors.destroy();
    vboColors.create();
    if(vboColors.isCreated())
        vboColors.bind();
    else
        return false;
    vboColors.setUsagePattern(QOpenGLBuffer::StaticDraw);
    program->enableAttributeArray(2);
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

	eboTriangles.destroy();
	eboTriangles.create();
	if(eboTriangles.isCreated())
		eboTriangles.bind();
	else
		return false;
	eboTriangles.setUsagePattern(QOpenGLBuffer::StaticDraw);

    fillVBOs(replicatedVertices, normals, perFaceTriangles, colors);

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

void TriangleMesh::render(QOpenGLFunctions &gl)
{
	vao.bind();
	eboTriangles.bind();
    gl.glDrawElements(GL_TRIANGLES, triangles.size(), GL_UNSIGNED_INT, nullptr);
	vao.release();
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

void TriangleMesh::fillVBOs(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles, vector<QVector3D> vertColors)
{
	vboVertices.bind();
	vboVertices.allocate(&replicatedVertices[0], 3 * sizeof(float) * replicatedVertices.size());
	vboVertices.release();

    vboColors.bind();
    vboColors.allocate(&vertColors[0], 3 * sizeof(float) * vertColors.size());
    vboColors.release();

	vboNormals.bind();
	vboNormals.allocate(&normals[0], 3 * sizeof(float) * normals.size());
	vboNormals.release();

	eboTriangles.bind();
	eboTriangles.allocate(&perFaceTriangles[0], sizeof(int) * perFaceTriangles.size());
	eboTriangles.release();
}


//Corner table calculation
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
        allCornersEdges.push_back(CornerEdge(e1Min, e1Max, i));
        allCornersEdges.push_back(CornerEdge(e2Min, e2Max, i+1));
        allCornersEdges.push_back(CornerEdge(e3Min, e3Max, i+2));

        //Assign Corner to Vertex
        cornerVertex[triangles[i*3]]=i;
        cornerVertex[triangles[i*3+1]]=i+1;
        cornerVertex[triangles[i*3+2]]=i+2;

        //Build Corner table
        cornersTable.resize(allCornersEdges.size());
        for (int i =0; i < allCornersEdges.size(); i++){
            CornerEdge e1 = allCornersEdges[i];
            for (int j =0; j < allCornersEdges.size(); j++){
                CornerEdge e2 = allCornersEdges[j];
                if (i!=j && CornerEdge::Compare(e1, e2)){
                    cornersTable[i] = e2;
                    cornersTable[j] = e1;
                }
            }

        }

    }

    std::cout<<"Corner table generated"<<std::endl;
}

vector<int> TriangleMesh::GetVertexNeighboors(int vert){
    vector<int> neighboors;
    //get corner of this vertex
    int vertCorner = cornerVertex[vert];
    int initialCorner = next(vertCorner);
    int nextCorner = initialCorner;
    //find vertices surrounding this corner
    while(true){
        int nextCorner = previous(nextCorner);
        if (nextCorner == -1) break; //non manifold mesh
        if (std::find(neighboors.begin(), neighboors.end(), nextCorner) == neighboors.end()) break;
        neighboors.push_back(nextCorner);
    }
    return neighboors;
}

vector<float> TriangleMesh::GaussianCurvature(float &min, float &max){
    vector<float> perVertCurvature(vertices.size());
    float maxCurvature = -9999999999;
    float minCurvature = 9999999999;

    //k = 2pi - sum(anglei) / area.
    for (int i = 0; i < vertices.size(); i++){
        //get neighboorhood of the vertex
        vector<int> neighboors = GetVertexNeighboors(i);
        //get sum of neighbooring angles
        float angle = 0.0;
        for (int j=0; j < neighboors.size(); j++){
            QVector3D v1;
            QVector3D v2;
            float a = acos(QVector3D::dotProduct(v1, v2));
            angle += a;
        }

        float curvature = 2*PI - angle;
        if (curvature > maxCurvature) maxCurvature = curvature;
        if (curvature < minCurvature) minCurvature = curvature;
        perVertCurvature[i] = curvature;
    }
    min = minCurvature; max = maxCurvature;
    return perVertCurvature;
}

void TriangleMesh::GetColors(vector<QVector3D> &vertColors, vector<float> vertCurvature, float min, float max){
    vertColors.resize(vertCurvature.size());
    for (unsigned int i = 0; i<vertColors.size(); i++){
        if (vertCurvature[i] > 0) vertColors[i] = QVector3D(vertCurvature[i]/max, 0, 0);
        else if (vertCurvature[i] < 0) vertColors[i] = QVector3D(vertCurvature[i]/-min, 0, 0);
        else vertColors[i] = QVector3D(0, 0, 0);
    }
}



