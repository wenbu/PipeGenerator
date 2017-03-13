#include "MeshUtil.h"
#include "ArrayUtil.h"
#include <maya/MIntArray.h>
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
