#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H


#include <vector>
#include <QVector3D>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#include "smoothingmatrix.h"


using namespace std;

struct CornerEdge
{
    int edgeMin, edgeMax, cornerVertex;
    CornerEdge(int eM, int eX, int v): edgeMin(eM), edgeMax(eX), cornerVertex(v){}
    CornerEdge(){}
    static bool Compare(CornerEdge e1, CornerEdge e2);
};

class TriangleMesh
{

public:
	TriangleMesh();
	void addVertex(const QVector3D &position);
	void addTriangle(int v0, int v1, int v2);

	void buildCube();

	bool init(QOpenGLShaderProgram *program);
	void destroy();

	void render(QOpenGLFunctions &gl);

    vector<int> cornersTable;
    vector<int> cornerVertex;

public:
    void DisplayGaussianCurvature();
    void DisplayMeanCurvature();
    void IteractiveSmoothing(int nSteps);
    void BiIteractiveSmoothing(int nSteps);
    void GlobalSmoothing(int percent);
    void DetMagnification(QVector3D l);

private:
	void buildReplicatedVertices(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles);
    vector<QVector3D> buildReplicatedColors(vector<QVector3D> currColors);
    void fillVBOs(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles);
    void updateColors(vector<QVector3D> &newColors);
    void updateVertices();

private:
    vector<QVector3D> vertices;
	vector<int> triangles;
    vector<QVector2D> textureCoordinates;
    pair<QVector3D, QVector3D> boundingBox;

    void buildCornerTable();
    vector<int> GetVertexNeighboors(unsigned int vert);

    vector<float> GaussianCurvature();
    vector<float> MeanCurvature();

    QVector3D ComputeLaplacian(int v_index, bool uniform);
    void IteractiveSmoothingStep();
    void BiIteractiveSmoothingStep();

    void buildSmoothingMatrix(int);
    void solveSparseSmoothing(int boundaryHeight);

    void getNoise(float l1, float l2, float l3);
    vector<QVector3D> SmoothingSteps(vector<QVector3D> &toSmooth, int n);

    void GetColors(vector<QVector3D> &vertColors, vector<float>&vertCurvatures);

    QOpenGLVertexArrayObject vao;
    QOpenGLBuffer vboVertices, vboNormals, vboColors, eboTriangles;

};


#endif // TRIANGLEMESH_H
