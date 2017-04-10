#pragma once

#include <maya/MFnMesh.h>
#include <maya/MPoint.h>
#include <maya/MTransformationMatrix.h>
#include <vector>

namespace MeshUtil
{
	MStatus transformPoints(MFnMesh& mesh, int* vertexIndices, int numVertices, MTransformationMatrix transformMatrix);
	MStatus extrudeWithTransform(MFnMesh& mesh, int faceIndex, int* vertexIndices, int numVertices, MTransformationMatrix transformMatrix);
	MStatus extrudeWithNewPositions(MFnMesh& mesh, int faceIndex, int* vertexIndices, int numVertices, std::vector<MPoint> newPositions);
	MStatus doExtrudeWithTranslate(MFnMesh& mesh, int& firstVertexIndex, int numVertices, int faceIndex, MPoint& scalePivot, MVector translateDirection, double translateMagnitude);
	MStatus doExtrudeWithUniformScale(MFnMesh& mesh, int& firstVertexIndex, int numVertices, int faceIndex, MPoint scalePivot, double scaleFactor);

	MStatus extrude(MFnMesh& mesh, int faceIndex);
	MStatus extrudeWithTransform(MFnMesh& mesh, int faceIndex, MTransformationMatrix transform);
	// A positive distance is parallel to the normal vector; a negative distance antiparallel.
	MStatus slideFace(const MObject& mesh, MFnMesh& fnMesh, int faceIndex, double distance);
	MStatus getCentroid(const MFnMesh& mesh, int faceIndex, MPoint* o_centroid);
}