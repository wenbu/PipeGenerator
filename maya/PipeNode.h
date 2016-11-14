#pragma once

#include <maya/MPxNode.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <vector>

class PipeNode : MPxNode
{
public:
	static const MString nodeName;
	static const MString nodeClassify;
	static const MTypeId nodeID;
	PipeNode();
	~PipeNode();

	MStatus compute(const MPlug& plug, MDataBlock& dataBlock);
	static void* creator();
	static MStatus initialize();
private:
	int numSides;
	double sectionRadius;
	double cornerAnglePerDivision;
	double pipeRotation;
	MString boltPath;
	double boltScale;
	int seed;
	bool createBeginningCap;
	bool createEndCap;
	MObject* boltShapeObjPtr;
	int numCapSides;
	double capThickness;
	double startCapTwist;
	double endCapTwist;
	double capRadius;
	double initialCapTightness;
	double secondCapTightness;
	double thirdCapTightness;
	double sideCapTightness;
	MFnMesh fnMesh;
	MFnMeshData dataCreator;
	MObject outputData;
	MObject objMesh;

	std::vector<MPoint> cvPoints;
	std::vector<MVector> spanVectors;
	std::vector<double> angles;
	std::vector<double> radii;
	std::vector<MPoint> centerPoints;
	std::vector<MVector> circleNormals;
	std::vector<double> spanLengthDeltas; // is this a good name? a "span length delta" is the distance from each CV to the point at which the arc intersects
	std::vector<int> cornerDivisions;
	std::vector<double> spanLengthBetweenCorners;
	std::vector<double> spanLengthFromLastCorner;

	MStatus getPipeData();
	MStatus makeBaseMesh();
	MStatus makePipe();
	MStatus makeCaps();
	int numVertsInPipe();
	MStatus transformPoints(int* vertexIndices, int numVertices, MTransformationMatrix transformMatrix);
	MStatus extrudeWithTransform(int faceIndex, int* vertexIndices, int numVertices, MTransformationMatrix transformMatrix);
	MStatus extrudeWithNewPositions(int faceIndex, int* vertexIndices, int numVertices, std::vector<MPoint> newPositions);
	MStatus doExtrudeWithTranslate(int& firstVertexIndex, int numVertices, int faceIndex, MPoint& scalePivot, MVector translateDirection, double translateMagnitude);
	MStatus doExtrudeWithUniformScale(int& firstVertexIndex, int numVertices, int faceIndex, MPoint scalePivot, double scaleFactor);

	// probably a better way to do this than pass in bools
	template<typename T> static MStatus registerNumericAttribute(MObject* attribute, const char* fullName, const char* shortName, MFnNumericData::Type type,
		T defaultValue = 0, bool hasMin = false, T minValue = 0, bool hasMax = false, T maxValue = 0);

	static MObject attrObjCurve;
	static MObject attrObjRadialSubdivisions;
	static MObject attrObjSectionRadius;
	static MObject attrObjRadiusArray;
	static MObject attrObjCornerAnglePerDivision;
	static MObject attrObjNumCapSides;
	static MObject attrObjCapThickness;
	static MObject attrObjBeginningCapRotation;
	static MObject attrObjEndCapRotation;
	static MObject attrObjCreateBeginningCap;
	static MObject attrObjCreateEndCap;
	static MObject attrObjCapRadius;
	static MObject attrObjInitialCapTightness;
	static MObject attrObjSecondCapTightness;
	static MObject attrObjThirdCapTightness;
	static MObject attrObjSideCapTightness;
	static MObject attrObjPipeRotation;
	static MObject attrObjBoltPath;
	static MObject attrObjBoltScale;
	static MObject attrObjSeed;
	static MObject attrObjMesh;
};

