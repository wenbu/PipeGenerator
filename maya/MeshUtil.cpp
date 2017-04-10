#include "MeshUtil.h"
#include "ArrayUtil.h"
#include <maya/MIntArray.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MMatrix.h>

MStatus MeshUtil::transformPoints(MFnMesh& fnMesh, int* vertexIndices, int numVertices, MTransformationMatrix transformMatrix)
{
	for (int i = 0; i < numVertices; i++)
	{
		int vidx = vertexIndices[i];
		MPoint pt;
		fnMesh.getPoint(vidx, pt, MSpace::kObject);
		fnMesh.setPoint(vidx, pt * transformMatrix.asMatrix(), MSpace::kObject);
		// error handling?
	}
	return MS::kSuccess;
}

MStatus MeshUtil::doExtrudeWithTranslate(MFnMesh& fnMesh, int& firstVertexIndex, int numVertices, int faceIndex, MPoint& scalePivot, MVector translateDirection, double translateMagnitude)
{
	int* vertexIndices = ArrayUtil::range(firstVertexIndex, firstVertexIndex + numVertices);
	MTransformationMatrix transform;
	MVector translateVector = translateDirection * translateMagnitude;
	transform.setTranslation(translateVector, MSpace::kWorld);
	extrudeWithTransform(fnMesh, faceIndex, vertexIndices, numVertices, transform);
	delete[] vertexIndices;

	firstVertexIndex += numVertices;
	scalePivot += translateVector;
	return MS::kSuccess;
}

MStatus MeshUtil::doExtrudeWithUniformScale(MFnMesh& fnMesh, int& firstVertexIndex, int numVertices, int faceIndex, MPoint scalePivot, double scaleFactor)
{
	int* vertexIndices = ArrayUtil::range(firstVertexIndex, firstVertexIndex + numVertices);
	MTransformationMatrix transform;
	double scale[3] = { scaleFactor, scaleFactor, scaleFactor };
	transform.setScale(scale, MSpace::kWorld);
	transform.setScalePivot(scalePivot, MSpace::kWorld, false);
	extrudeWithTransform(fnMesh, faceIndex, vertexIndices, numVertices, transform);
	delete[] vertexIndices;

	firstVertexIndex += numVertices;
	return MS::kSuccess;
}

MStatus MeshUtil::extrudeWithTransform(MFnMesh& fnMesh, int faceIndex, int* vertexIndices, int numVertices, MTransformationMatrix transformMatrix)
{
	MIntArray faceArray(1, faceIndex);
	fnMesh.extrudeFaces(faceArray, 1, &MFloatVector(0, 0, 0), true);

	for (int i = 0; i < numVertices; i++)
	{
		int vidx = vertexIndices[i];
		MPoint pt;
		fnMesh.getPoint(vidx, pt, MSpace::kObject);
		fnMesh.setPoint(vidx, pt * transformMatrix.asMatrix(), MSpace::kObject);
		// error handling?
	}
	return MS::kSuccess;
}

MStatus MeshUtil::extrudeWithNewPositions(MFnMesh& fnMesh, int faceIndex, int* vertexIndices, int numVertices, std::vector<MPoint> newPositions)
{
	MIntArray faceArray(1, faceIndex);
	fnMesh.extrudeFaces(faceArray, 1, &MFloatVector(0, 0, 0), true);
	// assume numVertices == newPositions.size()
	for (int i = 0; i < numVertices; i++)
	{
		int vidx = vertexIndices[i];
		MPoint npt = newPositions[i];
		fnMesh.setPoint(vidx, npt, MSpace::kObject);
		// error handling?
	}
	return MS::kSuccess;
}

MStatus MeshUtil::extrude(MFnMesh& fnMesh, int faceIndex)
{
	MIntArray faceIndices(1, faceIndex);
	fnMesh.extrudeFaces(faceIndices, 1, &MFloatVector(0, 0, 0), true);
	return MS::kSuccess;
}

MStatus MeshUtil::extrudeWithTransform(MFnMesh& fnMesh, int faceIndex, MTransformationMatrix transform)
{
	MIntArray faceIndices(1, faceIndex);

	// after extrusion, extruded face retains the same index, but its original vertices stay put
	// the newly extruded face has new vertices
	fnMesh.extrudeFaces(faceIndices, 1, &MFloatVector(0, 0, 0), true);

	MIntArray newVertexIndices;
	fnMesh.getPolygonVertices(faceIndex, newVertexIndices);
	for (int i = 0; i < newVertexIndices.length(); i++)
	{
		int v_idx = newVertexIndices[i];
		MPoint pt;
		fnMesh.getPoint(v_idx, pt, MSpace::kObject);
		fnMesh.setPoint(v_idx, pt * transform.asMatrix(), MSpace::kObject);
		// error handling?
	}
	return MS::kSuccess;
}

MStatus MeshUtil::slideFace(const MObject& mesh, MFnMesh& fnMesh, int faceIndex, double distance)
{
	// TODO
	MIntArray vertexIndices;
	fnMesh.getPolygonVertices(faceIndex, vertexIndices);
	MVector normal;
	fnMesh.getPolygonNormal(faceIndex, normal);

	MItMeshPolygon it_faces(mesh);
	int prevIndex; // don't care about this
	it_faces.setIndex(faceIndex, prevIndex);

	MIntArray connectedEdges;
	it_faces.getConnectedEdges(connectedEdges);

	for (int i = 0; i < connectedEdges.length(); i++)
	{
		int2 edgeVertices;
		fnMesh.getEdgeVertices(connectedEdges[i], edgeVertices);

		// TODO: get the edge vector closest to pointing away from normal,
		// project normal onto it
	}
	return MS::kSuccess;
}

MStatus MeshUtil::getCentroid(const MFnMesh& fnMesh, int faceIndex, MPoint* o_centroid)
{
	MIntArray vertexIndices;
	fnMesh.getPolygonVertices(faceIndex, vertexIndices);

	double4 centroid;
	for (int i = 0; i < vertexIndices.length(); i++)
	{
		MPoint pt;
		fnMesh.getPoint(vertexIndices[i], pt, MSpace::kObject);

		double4 arr_pt;
		pt.get(arr_pt);

		for (int d = 0; d < 4; d++)
		{
			centroid[d] += arr_pt[d];
		}
	}

	for (int d = 0; d < 4; d++)
	{
		centroid[d] /= vertexIndices.length();
	}

	*o_centroid = MPoint(centroid);
	return MS::kSuccess;
}