#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H


#include <vector>
#include <QVector3D>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>


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

public:
	void addVertex(const QVector3D &position);
	void addTriangle(int v0, int v1, int v2);

	void buildCube();

	bool init(QOpenGLShaderProgram *program);
	void destroy();

	void render(QOpenGLFunctions &gl);

    vector<CornerEdge> cornersTable;
    vector<int> cornerVertex;

private:
	void buildReplicatedVertices(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles);
	void fillVBOs(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles);

private:
	vector<QVector3D> vertices;
	vector<int> triangles;

    void buildCornerTable();
    vector<int> GetVertexNeighboors(int vert);
    void GaussianCurvature();

    QOpenGLVertexArrayObject vao;
	QOpenGLBuffer vboVertices, vboNormals, eboTriangles;

};


#endif // TRIANGLEMESH_H
