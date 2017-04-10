#pragma once

#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MPxNode.h>

class CapNode : MPxNode
{
public:
	static const MString nodeName;
	static const MString nodeClassify;
	static const MTypeId nodeID;
	CapNode();
	~CapNode();

	MStatus compute(const MPlug& plug, MDataBlock& dataBlock);
	static void* creator();
	static MStatus initialize();
private:
	MStatus makeCaps();

	// Return mean of distances from centroid to each vertex of specified face.
	// Placeholder until we do proper polygon skeleton stuff.
	MStatus getPseudoRadius(const MFnMesh& fnMesh, const MPoint& centroid, int faceIndex, double* o_pseudoRadius);

	MObject* boltShapeObjPtr;      // pointer to bolt object, once path is resolved
	MFnMeshData dataCreator;
	MObject outputData;
	MObject objMesh;

	MObject inputMesh;
	MString boltPath;              // path to bolt object
	double boltScale;              // scaling factor for bolts
	int seed;                      // PRNG seed
	int beginningCapIndex;         // index of beginning cap
	int endCapIndex;               // index of end cap; -1 if none
	int numCapSides;               // number of sides in cap; 0 = same as pipe
	double capThickness;           // thickness of caps
	double startCapTwist;          // twist of beginning cap
	double endCapTwist;            // twist of end cap
	double capRadius;              // radius of cap
	double initialCapTightness;    // edge loop tightness at cap-pipe joint
	double secondCapTightness;     // edge loop tightness at next joint
	double thirdCapTightness;      // edge loop tightness at next joint
	double sideCapTightness;       // edge loop tightness of sides

	static MObject attrObjMeshIn;
	static MObject attrObjNumCapSides;
	static MObject attrObjCapThickness;
	static MObject attrObjBeginningCapRotation;
	static MObject attrObjEndCapRotation;
	static MObject attrObjBeginningCapIndex;
	static MObject attrObjEndCapIndex;
	static MObject attrObjCapRadius;
	static MObject attrObjInitialCapTightness;
	static MObject attrObjSecondCapTightness;
	static MObject attrObjThirdCapTightness;
	static MObject attrObjSideCapTightness;
	static MObject attrObjBoltPath;
	static MObject attrObjBoltScale;
	static MObject attrObjSeed;
	static MObject attrObjMeshOut;
};