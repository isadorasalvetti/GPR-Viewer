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
    int vertexA, vertexB, corner;

    bool checkDouble(CornerEdge b){
        return (b.vertexA == vertexB) && (b.vertexB == vertexA);
    }

    bool operator<(const CornerEdge &cEdge) { return (vertexA < cEdge.vertexA) || ((vertexA == cEdge.vertexA) && (vertexB < cEdge.vertexB)); }
    bool operator==(const CornerEdge &cEdge) { return (vertexA == cEdge.vertexA) && (vertexB == cEdge.vertexB); }
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

private:
	void buildReplicatedVertices(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles);
	void fillVBOs(vector<QVector3D> &replicatedVertices, vector<QVector3D> &normals, vector<unsigned int> &perFaceTriangles);

private:
	vector<QVector3D> vertices;
	vector<int> triangles;

    vector<CornerEdge> tableCornerEdges;
    vector<int> tableCornerVertices;

    void buildCornerTable();
    vector<int> GetVertexNeighboors();
    void GaussianCurvature();

    QOpenGLVertexArrayObject vao;
	QOpenGLBuffer vboVertices, vboNormals, eboTriangles;

};


#endif // TRIANGLEMESH_H
