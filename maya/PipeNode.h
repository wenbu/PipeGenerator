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

