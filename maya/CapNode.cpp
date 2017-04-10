#include "CapNode.h"

#include "ArrayUtil.h"
#include "AttributeUtil.h"
#include "MeshUtil.h"
#include "TransformUtil.h"

#include <maya/MDagPath.h>
#include <maya/MFnStringData.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MMatrix.h>
#include <maya/MQuaternion.h>
#include <maya/MSelectionList.h>

#include <algorithm>

const MTypeId CapNode::nodeID(0xBEF0);
const MString CapNode::nodeName("capnode");
const MString CapNode::nodeClassify("utility/general");

MObject CapNode::attrObjNumCapSides;
MObject CapNode::attrObjCapThickness;
MObject CapNode::attrObjBeginningCapRotation;
MObject CapNode::attrObjEndCapRotation;
MObject CapNode::attrObjBeginningCapIndex;
MObject CapNode::attrObjEndCapIndex;
MObject CapNode::attrObjCapRadius;
MObject CapNode::attrObjInitialCapTightness;
MObject CapNode::attrObjSecondCapTightness;
MObject CapNode::attrObjThirdCapTightness;
MObject CapNode::attrObjSideCapTightness;
MObject CapNode::attrObjBoltPath;
MObject CapNode::attrObjBoltScale;
MObject CapNode::attrObjSeed;

CapNode::CapNode()
{
}

CapNode::~CapNode()
{
}

MStatus CapNode::compute(const MPlug& plug, MDataBlock& dataBlock)
{
	MStatus stat;
	MDataHandle dataHandleInMesh = dataBlock.inputValue(attrObjMeshIn, &stat); if (!stat) return stat;

	MDataHandle dataHandleNumCapSides = dataBlock.inputValue(attrObjNumCapSides, &stat); if (!stat) return stat;
	MDataHandle dataHandleCapThickness = dataBlock.inputValue(attrObjCapThickness, &stat); if (!stat) return stat;
	MDataHandle dataHandleBeginningCapRotation = dataBlock.inputValue(attrObjBeginningCapRotation, &stat); if (!stat) return stat;
	MDataHandle dataHandleEndCapRotation = dataBlock.inputValue(attrObjEndCapRotation, &stat); if (!stat) return stat;
	MDataHandle dataHandleCapRadius = dataBlock.inputValue(attrObjCapRadius, &stat); if (!stat) return stat;

	MDataHandle dataHandleCreateBeginningCap = dataBlock.inputValue(attrObjBeginningCapIndex, &stat); if (!stat) return stat;
	MDataHandle dataHandleCreateEndCap = dataBlock.inputValue(attrObjEndCapIndex, &stat); if (!stat) return stat;

	MDataHandle dataHandleInitialCapTightness = dataBlock.inputValue(attrObjInitialCapTightness, &stat); if (!stat) return stat;
	MDataHandle dataHandleSecondCapTightness = dataBlock.inputValue(attrObjSecondCapTightness, &stat); if (!stat) return stat;
	MDataHandle dataHandleThirdCapTightness = dataBlock.inputValue(attrObjThirdCapTightness, &stat); if (!stat) return stat;
	MDataHandle dataHandleSideCapTightness = dataBlock.inputValue(attrObjSideCapTightness, &stat); if (!stat) return stat;

	MDataHandle dataHandleBoltPath = dataBlock.inputValue(attrObjBoltPath, &stat); if (!stat) return stat;
	MDataHandle dataHandleBoltScale = dataBlock.inputValue(attrObjBoltScale, &stat); if (!stat) return stat;
	MDataHandle dataHandleSeed = dataBlock.inputValue(attrObjSeed, &stat); if (!stat) return stat;

	MDataHandle dataHandleOutMesh = dataBlock.outputValue(attrObjMeshOut, &stat); if (!stat) return stat;

	// asMesh or asMeshTransformed?
	inputMesh = dataHandleInMesh.asMeshTransformed();

	boltPath = dataHandleBoltPath.asString();
	boltScale = dataHandleBoltScale.asDouble();
	seed = dataHandleSeed.asInt();
	beginningCapIndex = dataHandleCreateBeginningCap.asInt();
	endCapIndex = dataHandleCreateEndCap.asInt();
	srand(seed);

	if (boltPath.length() == 0)
	{
		boltShapeObjPtr = NULL;
	}
	else
	{
		MSelectionList boltSelectionList;
		MDagPath dagPath;
		stat = boltSelectionList.add(boltPath);
		if (!stat) return stat;
		stat = boltSelectionList.getDagPath(0, dagPath);
		if (!stat) return stat;
		MObject node = dagPath.node(&stat);
		if (!stat) return stat;
		if (!node.hasFn(MFn::kMesh))
		{
			// if transform, try extending to shape and see if it's a mesh
			stat = dagPath.extendToShape();
			if (!stat) return stat;
			node = dagPath.node(&stat);
			if (!stat) return stat;
			if (!node.hasFn(MFn::kMesh))
			{
				// still not a mesh; give up
				boltShapeObjPtr = NULL;
			}
		}
		boltShapeObjPtr = &dagPath.node(&stat);
		if (!stat) return stat;
	}

	numCapSides = dataHandleNumCapSides.asInt();
	capThickness = dataHandleCapThickness.asDouble();
	startCapTwist = dataHandleBeginningCapRotation.asDouble();
	endCapTwist = dataHandleEndCapRotation.asDouble();
	capRadius = dataHandleCapRadius.asDouble();

	initialCapTightness = dataHandleInitialCapTightness.asDouble();
	secondCapTightness = dataHandleSecondCapTightness.asDouble();
	thirdCapTightness = dataHandleThirdCapTightness.asDouble();
	sideCapTightness = dataHandleSideCapTightness.asDouble();

	// force triangular cap if numCapSides < 3
	// numCapSides = 0 is a special case that indicates cylinder so 
	// don't force for numCapSides = 0
	if (numCapSides < 0 || numCapSides == 1 || numCapSides == 2)
	{
		numCapSides = 3;
	}

	return MStatus::kSuccess;
}

// Input is a mesh and two face indices. Create caps on these faces.
MStatus CapNode::makeCaps()
{
	int faceIndices[2] = { beginningCapIndex, endCapIndex }; // initial indices of target faces
	//int faces[2] = { 1, 0 };                                 // ??? what the fuck is this?
	//int dataIndices[2] = { 0, -1 };                          // ??? what the fuck is this?
	//int extrusionDirs[2] = { -1, 1 };                        // seems hacky; consider refactor
	bool shouldFlipVerts[2] = { true, false };               // hack to get around weird maya behavior; should see if still needed
	//int pipeVertexOffsets[2] = { 0, numVertsInPipe() - numSides }; // vertex index offsets for target faces; should consider refactor
	double twists[2] = { startCapTwist, endCapTwist };       // twist value, in degrees

	//int vertexOffset = numVertsInPipe();

	// overall flow:
	// (note: 'down' means 'in direction of face normal'; 'up' is antiparallel to this
	//
	// for each target face:
	// 1) extrude
	// 2) if numCapSides == 0
	//      scale outwards; would like to do something with polygon skeleton but for now uniform scale suffices
	//    else
	//      move verts appropriately & scale outwards according to cap size
	// 3) extrude & move down according to cap thickness
	// 4) cap off
	//
	// edge loops as needed (i.e. before each step)

	// currently hardcoded to two faces; can probably generalize to n faces later
	for (int i = 0; i < 2; i++)
	{
		int capIndex = faceIndices[i];
		if (capIndex == -1) continue;

		MFnMesh fnMesh(inputMesh);
		MPoint scalePivot;
		fnMesh.getPoint(capIndex, scalePivot);
		bool shouldFlipVert = shouldFlipVerts[i];
		//int pipeVertexOffset = pipeVertexOffsets[i];
		//int dataIndex = dataIndices[i];
		//int extrusionDir = extrusionDirs[i];
		//int face = faces[i];
		double twist = twists[i];
		MVector normal;
		fnMesh.getPolygonNormal(capIndex, normal); // needs normalization?

		/*
		 * need to move polygon back (i.e. opposite direction of normal) by a distance = cap thickness
		 * don't just slide back along normal; vertices should be slid back along their edges such
		 * that the polygon normal remains the same (polygon area may change if these incident edges
		 * are not parallel to polygon normal)
		 */

		// BACK UP FOR EDGE LOOP ONE HERE
		MeshUtil::slideFace(objMesh, fnMesh, capIndex, -initialCapTightness);
		/*
		// back up existing verts a bit for edge loop
		int* vertexIndices = ArrayUtil::range(pipeVertexOffset, pipeVertexOffset + numSides);
		MTransformationMatrix transform;
		double capThick = shouldCreateCap ? capThickness : 0;
		MVector translateVector = normal * (-extrusionDir * (capThick + initialCapTightness));
		transform.setTranslation(translateVector, MSpace::kWorld);
		MeshUtil::transformPoints(fnMesh, vertexIndices, numSides, transform);
		delete[] vertexIndices;
		scalePivot += translateVector;
		*/

		// extrude down for initial edge loop 1
		//MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, normal, extrusionDir * initialCapTightness);
		MeshUtil::extrude(fnMesh, capIndex);
		MeshUtil::slideFace(fnMesh, capIndex, initialCapTightness);

		//if (!shouldCreateCap) continue;

		MPoint centroid;
		MeshUtil::getCentroid(fnMesh, capIndex, &centroid);

		double pseudoRadius;
		getPseudoRadius(fnMesh, centroid, capIndex, &pseudoRadius);
		double extrudeAmount = (pseudoRadius + initialCapTightness) / pseudoRadius;

		// extrude out for initial edge loop 2
		//MeshUtil::doExtrudeWithUniformScale(fnMesh, vertexOffset, numSides, face, scalePivot, (sectionRadius + initialCapTightness) / sectionRadius);
		MeshUtil::extrudeWithTransform(fnMesh, capIndex, TransformUtil::getUniformScale(pseudoRadius, centroid));

		// make cap shape; 0 = circle so do nothing
		if (numCapSides > 0)
		{
			// consider pulling this logic into a separate method

			int* vertexIndices = ArrayUtil::range(vertexOffset, vertexOffset + numSides);

			// get around weird maya behavior
			if (shouldFlipVert)
			{
				std::swap(vertexIndices[0], vertexIndices[1]);
			}

			MTransformationMatrix capTransform;
			MVector vRotSrc(0, 0, 1);
			MVector vRotDst = normal * extrusionDir;
			MQuaternion quat(vRotSrc, vRotDst);
			capTransform.rotateTo(quat);
			MPoint vTrnSrc(0, 0, 0);
			MPoint vTrnDst(scalePivot);
			MVector vTrn = vTrnDst - vTrnSrc;
			capTransform.setTranslation(vTrn, MSpace::kWorld);
			double scaleFactor = capRadius - secondCapTightness;
			double scale[3] = { scaleFactor, scaleFactor, scaleFactor };
			capTransform.setScale(scale, MSpace::kWorld);

			MMatrix transformMatrix = capTransform.asMatrix();

			// build reference polygon at origin
			std::vector<MPoint> refPoints;
			std::vector<MPoint> boltPoints;
			double dtheta = 2 * M_PI / numCapSides;
			for (int s = 0; s < numCapSides; s++)
			{
				double x = cos(s * dtheta + twist * (2 * M_PI / 360));
				double y = sin(s * dtheta + twist * (2 * M_PI / 360));
				double z = 0;
				refPoints.push_back(MPoint(x, y, z) * transformMatrix);

				double boltScale = (capRadius + sectionRadius) / (2 * capRadius);
				boltPoints.push_back(MPoint(boltScale*x, boltScale*y, z) * transformMatrix);
			}

			std::vector<MVector> polyVectors;
			std::vector<MPoint> rotatedRefPoints(numCapSides);
			std::rotate_copy(refPoints.begin(), refPoints.begin() + 1, refPoints.end(), rotatedRefPoints.begin());
			for (int j = 0; j < numCapSides; j++)
			{
				polyVectors.push_back(rotatedRefPoints[j] - refPoints[j]);
			}

			// figure out where cap verts go
			int numPtsRemaining = numSides;
			int numPtsOnThisSide = 0;
			std::vector<MPoint> newPositions;
			newPositions.reserve(numSides);
			for (int vidx = 0; vidx < numCapSides; vidx++)
			{
				MVector v = polyVectors[vidx];
				numPtsOnThisSide = numPtsRemaining / (numCapSides - vidx);
				numPtsRemaining -= numPtsOnThisSide;
				std::vector<MPoint> newPositionsOnThisSide;
				newPositionsOnThisSide.resize(numPtsOnThisSide);
				double edgeScaleFactor = sideCapTightness / v.length();
				for (int pidx = 0; pidx < numPtsOnThisSide; pidx++)
				{
					// pidx == 0 -> vert on polygon vertex
					// pidx == 1 -> vert on near edge loop
					// pidx == 2 -> vert on far edge loop
					// all others -> equally spaced between near and far edge loops
					double f = 0;
					int nidx = 0;
					if (pidx == 0)
					{
						// do nothing
					}
					else if (pidx == 1)
					{
						f = edgeScaleFactor;
						nidx = 1;
					}
					else if (pidx == 2)
					{
						f = 1 - edgeScaleFactor;
						nidx = numPtsOnThisSide - 1;
					}
					else
					{
						int numExcessPtsOnThisSide = numPtsOnThisSide - 3;
						int eidx = pidx - 3;
						double numerator = eidx + 1.0;
						double denominator = numExcessPtsOnThisSide + 1;
						double midSectionLength = v.length() - 2 * sideCapTightness;
						double midSectionFactor = midSectionLength / v.length();
						f = edgeScaleFactor + (numerator / denominator) * midSectionFactor;
						nidx = pidx - 1;
					}
					newPositionsOnThisSide[nidx] = refPoints[vidx] + (v * f);
				}
				newPositions.insert(newPositions.end(), newPositionsOnThisSide.begin(), newPositionsOnThisSide.end());
			}

			// rotate vertexIndices to minimize distortion from rotation
			double shortestDistance = std::numeric_limits<double>::max();
			int amountToRotate = 0;
			for (int j = 0; j < numSides; j++)
			{
				int* testIndices = ArrayUtil::rotateArrayCopy(vertexIndices, numSides, j);
				for (int k = 0; k < numSides; k++)
				{
					testIndices[k] = testIndices[k] - numSides;
				}
				double distance = 0;

				for (int iidx = 0; iidx < numSides; iidx++)
				{
					int pidx = testIndices[iidx];
					MPoint meshPt;
					fnMesh.getPoint(pidx, meshPt, MSpace::kObject);
					MPoint capPt = newPositions[iidx];
					distance += (capPt - meshPt).length();
				}
				if (distance < shortestDistance)
				{
					shortestDistance = distance;
					amountToRotate = j;
				}
			}
			if (amountToRotate != 0)
			{
				int* rotatedVertexIndices = ArrayUtil::rotateArrayCopy(vertexIndices, numSides, amountToRotate);
				delete[] vertexIndices;
				vertexIndices = rotatedVertexIndices;
			}
			MeshUtil::extrudeWithNewPositions(fnMesh, face, vertexIndices, numSides, newPositions);
			delete[] vertexIndices;
			vertexOffset += numSides;
		}

		// extrude out for second edge loop 1
		MeshUtil::doExtrudeWithUniformScale(fnMesh, vertexOffset, numSides, face, scalePivot, capRadius / (capRadius - secondCapTightness));

		// extrude down for second edge loop 2
		MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, normal, extrusionDir * secondCapTightness);

		// body extrusion
		MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, normal, extrusionDir * (capThickness - secondCapTightness - thirdCapTightness));

		// extrude down for third edge loop 1
		MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, normal, extrusionDir * thirdCapTightness);

		// extrude in for third edge loop 2
		MeshUtil::doExtrudeWithUniformScale(fnMesh, vertexOffset, numSides, face, scalePivot, (capRadius - thirdCapTightness) / capRadius);

		// extrude in
		MeshUtil::doExtrudeWithUniformScale(fnMesh, vertexOffset, numSides, face, scalePivot, 0.5);

		// bolts TODO
	}
	return MS::kSuccess;
}

MStatus CapNode::getPseudoRadius(const MFnMesh& fnMesh, const MPoint& centroid, int faceIndex, double* o_pseudoRadius)
{
	MIntArray vertexIndices;
	fnMesh.getPolygonVertices(faceIndex, vertexIndices);

	double meanDistance = 0;
	for (int i = 0; i < vertexIndices.length(); i++)
	{
		MPoint pt;
		fnMesh.getPoint(vertexIndices[i], pt, MSpace::kObject);

		meanDistance += pt.distanceTo(centroid);
	}

	meanDistance /= vertexIndices.length();
	*o_pseudoRadius = meanDistance;
	return MS::kSuccess;
}

void* CapNode::creator()
{
	return new CapNode;
}

MStatus CapNode::initialize()
{
	MFnTypedAttribute typedAttr;

	attrObjMeshOut = typedAttr.create("meshOut", "mo", MFnData::kMesh);
	typedAttr.setKeyable(false);
	typedAttr.setWritable(false);
	addAttribute(attrObjMeshOut);

	attrObjMeshIn = typedAttr.create("meshIn", "mi", MFnData::kMesh);
	typedAttr.setKeyable(false);
	typedAttr.setReadable(false);
	addAttribute(attrObjMeshIn);
	attributeAffects(attrObjMeshIn, attrObjMeshOut);

	AttributeUtil::registerNumericAttribute<int>(&attrObjNumCapSides, &attrObjMeshOut, "capShape", "cs", MFnNumericData::kInt, 4);
	AttributeUtil::registerNumericAttribute<double>(&attrObjCapThickness, &attrObjMeshOut, "capThickness", "ct", MFnNumericData::kDouble, 1.0);
	AttributeUtil::registerNumericAttribute<int>(&attrObjBeginningCapIndex, &attrObjMeshOut, "beginningCapIndex", "ibc", MFnNumericData::kBoolean, true);
	AttributeUtil::registerNumericAttribute<double>(&attrObjBeginningCapRotation, &attrObjMeshOut, "beginningCapRotation", "bcr", MFnNumericData::kDouble, 0.0, true, 0.0, true, 360.0);
	AttributeUtil::registerNumericAttribute<int>(&attrObjEndCapIndex, &attrObjMeshOut, "endCapIndex", "iec", MFnNumericData::kBoolean, true);
	AttributeUtil::registerNumericAttribute<double>(&attrObjEndCapRotation, &attrObjMeshOut, "endCapRotation", "ecr", MFnNumericData::kDouble, 0.0, true, 0.0, true, 360.0);
	AttributeUtil::registerNumericAttribute<double>(&attrObjCapRadius, &attrObjMeshOut, "capRadius", "cra", MFnNumericData::kDouble, 2.0, true, 0.0);
	AttributeUtil::registerNumericAttribute<double>(&attrObjInitialCapTightness, &attrObjMeshOut, "initialCapTightness", "ict", MFnNumericData::kDouble, 0.1, true, 0.001);
	AttributeUtil::registerNumericAttribute<double>(&attrObjSecondCapTightness, &attrObjMeshOut, "secondCapTightness", "sct", MFnNumericData::kDouble, 0.1, true, 0.001);
	AttributeUtil::registerNumericAttribute<double>(&attrObjThirdCapTightness, &attrObjMeshOut, "thirdCapTightness", "tct", MFnNumericData::kDouble, 0.1, true, 0.001);
	AttributeUtil::registerNumericAttribute<double>(&attrObjSideCapTightness, &attrObjMeshOut, "sideCapTightness", "sict", MFnNumericData::kDouble, 0.1, true, 0.001);
	attrObjBoltPath = typedAttr.create("boltPath", "bp", MFnData::kString, MFnStringData().create());
	typedAttr.setKeyable(true);
	typedAttr.setStorable(true);
	addAttribute(attrObjBoltPath);
	attributeAffects(attrObjBoltPath, attrObjMeshOut);

	AttributeUtil::registerNumericAttribute<double>(&attrObjBoltScale, &attrObjMeshOut, "boltScale", "bs", MFnNumericData::kDouble, 1.0);
	AttributeUtil::registerNumericAttribute<int>(&attrObjSeed, &attrObjMeshOut, "seed", "s", MFnNumericData::kInt, 2938);
}