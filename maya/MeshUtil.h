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
}