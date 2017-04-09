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
	MFnMesh fnMesh;
	MFnMeshData dataCreator;
	MObject outputData;
	MObject objMesh;

	// ---------- populated directly from attributes ----------
	int numSides;                  // number of sides of polygonal cross-section
	double sectionRadius;          // radius of cross-section
	double cornerAnglePerDivision; // angle subtended by each corner segment
	double pipeRotation;           // twist

	// ---------- derived from attributes ----------
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

	// ---------- functions ----------
	MStatus getPipeData();
	MStatus makeBaseMesh();
	MStatus makePipe();
	int numVertsInPipe();

	// ---------- attributes ----------
	static MObject attrObjCurve;
	static MObject attrObjRadialSubdivisions;
	static MObject attrObjSectionRadius;
	static MObject attrObjRadiusArray;
	static MObject attrObjCornerAnglePerDivision;
	static MObject attrObjPipeRotation;
	static MObject attrObjMesh;
};

