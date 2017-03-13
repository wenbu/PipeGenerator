#include "PipeNode.h"
#include "ArrayUtil.h"
#include "AttributeUtil.h"
#include "MeshUtil.h"
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MItCurveCV.h>
#include <maya/MPoint.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MQuaternion.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MFnStringData.h>
#include <algorithm>

const MTypeId PipeNode::nodeID(0xBEEF);
const MString PipeNode::nodeName("pipenode");
const MString PipeNode::nodeClassify("utility/general");

MObject PipeNode::attrObjCurve;
MObject PipeNode::attrObjRadialSubdivisions;
MObject PipeNode::attrObjSectionRadius;
MObject PipeNode::attrObjRadiusArray;
MObject PipeNode::attrObjCornerAnglePerDivision;
MObject PipeNode::attrObjNumCapSides;
MObject PipeNode::attrObjCapThickness;
MObject PipeNode::attrObjBeginningCapRotation;
MObject PipeNode::attrObjEndCapRotation;
MObject PipeNode::attrObjCreateBeginningCap;
MObject PipeNode::attrObjCreateEndCap;
MObject PipeNode::attrObjCapRadius;
MObject PipeNode::attrObjInitialCapTightness;
MObject PipeNode::attrObjSecondCapTightness;
MObject PipeNode::attrObjThirdCapTightness;
MObject PipeNode::attrObjSideCapTightness;
MObject PipeNode::attrObjPipeRotation;
MObject PipeNode::attrObjBoltPath;
MObject PipeNode::attrObjBoltScale;
MObject PipeNode::attrObjSeed;
MObject PipeNode::attrObjMesh;

PipeNode::PipeNode()
{
}

PipeNode::~PipeNode()
{
}

MStatus PipeNode::compute(const MPlug& plug, MDataBlock& dataBlock)
{
	MStatus stat;
	MDataHandle dataHandleInCurve = dataBlock.inputValue(attrObjCurve, &stat); if (!stat) return stat;
	MDataHandle dataHandleInRadialSubdivisions = dataBlock.inputValue(attrObjRadialSubdivisions, &stat); if (!stat) return stat;
	MDataHandle dataHandleSectionRadius = dataBlock.inputValue(attrObjSectionRadius, &stat); if (!stat) return stat;
	MDataHandle dataHandleCornerAnglePerDivision = dataBlock.inputValue(attrObjCornerAnglePerDivision, &stat); if (!stat) return stat;

	MDataHandle dataHandleNumCapSides = dataBlock.inputValue(attrObjNumCapSides, &stat); if (!stat) return stat;
	MDataHandle dataHandleCapThickness = dataBlock.inputValue(attrObjCapThickness, &stat); if (!stat) return stat;
	MDataHandle dataHandleBeginningCapRotation = dataBlock.inputValue(attrObjBeginningCapRotation, &stat); if (!stat) return stat;
	MDataHandle dataHandleEndCapRotation = dataBlock.inputValue(attrObjEndCapRotation, &stat); if (!stat) return stat;
	MDataHandle dataHandleCapRadius = dataBlock.inputValue(attrObjCapRadius, &stat); if (!stat) return stat;

	MDataHandle dataHandleCreateBeginningCap = dataBlock.inputValue(attrObjCreateBeginningCap, &stat); if (!stat) return stat;
	MDataHandle dataHandleCreateEndCap = dataBlock.inputValue(attrObjCreateEndCap, &stat); if (!stat) return stat;

	MDataHandle dataHandleInitialCapTightness = dataBlock.inputValue(attrObjInitialCapTightness, &stat); if (!stat) return stat;
	MDataHandle dataHandleSecondCapTightness = dataBlock.inputValue(attrObjSecondCapTightness, &stat); if (!stat) return stat;
	MDataHandle dataHandleThirdCapTightness = dataBlock.inputValue(attrObjThirdCapTightness, &stat); if (!stat) return stat;
	MDataHandle dataHandleSideCapTightness = dataBlock.inputValue(attrObjSideCapTightness, &stat); if (!stat) return stat;

	MDataHandle dataHandlePipeRotation = dataBlock.inputValue(attrObjPipeRotation, &stat); if (!stat) return stat;

	MDataHandle dataHandleBoltPath = dataBlock.inputValue(attrObjBoltPath, &stat); if (!stat) return stat;
	MDataHandle dataHandleBoltScale = dataBlock.inputValue(attrObjBoltScale, &stat); if (!stat) return stat;
	MDataHandle dataHandleSeed = dataBlock.inputValue(attrObjSeed, &stat); if (!stat) return stat;

	MDataHandle dataHandleOutMesh = dataBlock.outputValue(attrObjMesh, &stat); if (!stat) return stat;

	MObject curve = dataHandleInCurve.asNurbsCurveTransformed();
	numSides = dataHandleInRadialSubdivisions.asInt();
	sectionRadius = dataHandleSectionRadius.asDouble();
	cornerAnglePerDivision = dataHandleCornerAnglePerDivision.asDouble() * (2 * M_PI / 360);
	pipeRotation = dataHandlePipeRotation.asDouble() * (2 * M_PI / 360);
	boltPath = dataHandleBoltPath.asString();
	boltScale = dataHandleBoltScale.asDouble();
	seed = dataHandleSeed.asInt();
	createBeginningCap = dataHandleCreateBeginningCap.asBool();
	createEndCap = dataHandleCreateEndCap.asBool();

	srand(seed);
	
	cvPoints.clear();
	spanVectors.clear();
	angles.clear();
	radii.clear();
	centerPoints.clear();
	circleNormals.clear();
	spanLengthDeltas.clear();
	cornerDivisions.clear();
	spanLengthBetweenCorners.clear();
	spanLengthFromLastCorner.clear();

	MItCurveCV cvIterator(curve);
	while (!cvIterator.isDone())
	{
		cvPoints.push_back(cvIterator.position(MSpace::kObject, &stat));
		if (!stat) return stat;
		stat = cvIterator.next();
		if (!stat) return stat;
	}

	MPlug radiusArrayPlug(thisMObject(), attrObjRadiusArray);
	for (int i = 0; i < cvPoints.size() - 2; i++)
	{
		stat = radiusArrayPlug.selectAncestorLogicalIndex(i, attrObjRadiusArray);
		if (!stat) return stat;
		radii.push_back(radiusArrayPlug.asFloat());
	}

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

	outputData = dataCreator.create(&stat);
	if (!stat) return stat;

	stat = getPipeData();
	if (!stat) return stat;
	stat = makeBaseMesh();
	if (!stat) return stat;
	stat = makePipe();
	if (!stat) return stat;
	stat = makeCaps();
	if (!stat) return stat;

	stat = fnMesh.updateSurface();
	if (!stat) return stat;
	stat = dataHandleOutMesh.setMObject(outputData);
	if (!stat) return stat;
	dataHandleOutMesh.setClean();

	return MStatus::kSuccess;
}

MStatus PipeNode::getPipeData()
{
	for (int i = 0; i < cvPoints.size() - 1; i++)
	{
		MVector spanVector(cvPoints[i + 1] - cvPoints[i]);
		spanVectors.push_back(spanVector);
	}

	for (int i = 0; i < spanVectors.size() - 1; i++)
	{
		MVector v1 = spanVectors[i] * -1;
		MVector v2 = spanVectors[i + 1];
		angles.push_back(v1.angle(v2));
	}

	// compute arc data for each CV
	// each nonterminal CV will have an arc such that the ends of the arc
	// are tangent to the curve (curve is assumed degree 1) and the radius
	// is as specified in the array
	for (int i = 0; i < angles.size(); i++)
	{
		MVector spanVector1 = spanVectors[i].normal() * -1;
		MVector spanVector2 = spanVectors[i + 1].normal();
		MPoint vertexPoint = cvPoints[i + 1];

		// build bisector
		MVector p1 = vertexPoint + spanVector1;
		MVector p2 = vertexPoint + spanVector2;
		MVector pm = (p1 + p2) / 2;
		MVector bisector = (pm - vertexPoint).normal();

		// distance from vert to arc center
		double theta = angles[i] / 2;
		double dist = radii[i] / sin(theta);

		// compute center point of arc
		centerPoints.push_back(vertexPoint + (bisector * dist));

		// get normal
		circleNormals.push_back(spanVector1 ^ spanVector2);

		// distance from the vertex to the point where the arc intersects
		spanLengthDeltas.push_back(radii[i] / tan(theta));

		// compute num segments
		// two metrics: angle-based and arclength based; take max
		// angle controlled by incoming attribute
		// arclength controlled by comparing segment length vs cross section side length
		// idea is to keep quads around the corner as square as possible overall
		double numDivisionsAngle = (M_PI - angles[i]) / cornerAnglePerDivision;

		double lengthOfCrossSectionSide = 2 * sectionRadius * sin(M_PI / numSides);
		double arcLength = (2 * M_PI * radii[i]) * (angles[i] / (2 * M_PI));
		double numDivisionsLength = arcLength / (8 * lengthOfCrossSectionSide);
		cornerDivisions.push_back(std::max(numDivisionsLength, numDivisionsAngle));
	}

	for (int i = 0; i < spanVectors.size(); i++)
	{
		double currentLength = spanVectors[i].length();

		if (i == 0)
		{
			// first item
			spanLengthBetweenCorners.push_back(currentLength - spanLengthDeltas[i]);
			spanLengthFromLastCorner.push_back(currentLength - spanLengthDeltas[i]);
		}
		else if (i == spanVectors.size() - 1)
		{
			// last item
			spanLengthBetweenCorners.push_back(currentLength - spanLengthDeltas[i - 1]);
			spanLengthFromLastCorner.push_back(currentLength);
		}
		else
		{
			// all other cases
			spanLengthBetweenCorners.push_back(currentLength - spanLengthDeltas[i] - spanLengthDeltas[i - 1]);
			spanLengthFromLastCorner.push_back(currentLength - spanLengthDeltas[i]);
		}
	}
	return MS::kSuccess;
}

MStatus PipeNode::makeBaseMesh()
{
	double radiansPerSide = (2 * M_PI) / numSides;

	std::vector<MPoint> untransformed;
	for (int i = 0; i < numSides; i++)
	{
		double x = sectionRadius * cos(pipeRotation + i * radiansPerSide);
		double y = sectionRadius * sin(pipeRotation + i * radiansPerSide);
		untransformed.push_back(MPoint(x, y, 0));
	}

	// build transform matrix to move points from origin to curve start
	// and orient to be perpendicular to curve
	MTransformationMatrix transform;
	MPoint targetPoint = cvPoints[0];
	MVector translateVector(targetPoint - MPoint(0, 0, 0));
	transform.setTranslation(translateVector, MSpace::kWorld);
	MVector targetNormal = spanVectors[0].normal();
	MQuaternion rotationQuat(MVector(0, 0, 1), targetNormal);
	transform.rotateTo(rotationQuat);
	MMatrix transformMatrix = transform.asMatrix();
	MPoint* points = new MPoint[numSides];
	for (int i = 0; i < numSides; i++)
	{
		points[i]  = (MPoint(untransformed[i] * transformMatrix));
	}

	int numPolygons = 1;
	MPointArray vertexArray(points, numSides);
	MIntArray polygonCounts(1, numSides); // initial size, initial value
	MIntArray polygonConnects(numSides);
	for (int i = 0; i < numSides; i++)
	{
		polygonConnects.set(i, i);
	}

	objMesh = fnMesh.create(numSides, numPolygons, vertexArray, polygonCounts, polygonConnects, outputData);

	delete[] points;
	return MS::kSuccess;
}

MStatus PipeNode::makePipe()
{
	MIntArray faceArray(1, 0);

	/*
	 * When we extrude a face, the extruded face maintains the same face index
	 * but the vertices are new, so we must update this offset every time we extrude
	 */
	int vertIndexStart = 0;
	for (int i = 0; i < spanVectors.size(); i++)
	{
		MVector span = spanVectors[i];
		MFloatVector spanVec(span.normal() * spanLengthBetweenCorners[i]);
		fnMesh.extrudeFaces(faceArray, 1, &spanVec, true);

		vertIndexStart += numSides;

		if (i < angles.size())
		{
			double theta = angles[i] - M_PI;
			double dtheta = theta / cornerDivisions[i];

			MTransformationMatrix transform;
			transform.setRotatePivot(centerPoints[i], MSpace::kWorld, false);
			transform.setToRotationAxis(circleNormals[i], dtheta);

			for (int j = 0; j < cornerDivisions[i]; j++)
			{
				// extrude but don't transform -- we will move points ourselves
				fnMesh.extrudeFaces(faceArray, 1, &MFloatVector(0, 0, 0), true);

				vertIndexStart += numSides;
				int* vertIndices = new int[numSides];
				for (int k = 0; k < numSides; k++)
				{
					vertIndices[k] = vertIndexStart + k;
				}
				MeshUtil::transformPoints(fnMesh, vertIndices, numSides, transform);
				delete[] vertIndices;
			}
		}
	}
	return MS::kSuccess;
}

MStatus PipeNode::makeCaps()
{
	// this is pretty gross
	bool shouldCreateCaps[2] = { createBeginningCap, createEndCap };
	int faces[2] = { 1, 0 };
	int dataIndices[2] = { 0, -1 };
	int extrusionDirs[2] = { -1, 1 };
	bool shouldFlipVerts[2] = { true, false };
	int pipeVertexOffsets[2] = { 0, numVertsInPipe() - numSides };
	double twists[2] = { startCapTwist, endCapTwist };
	
	int vertexOffset = numVertsInPipe();
	for (int i = 0; i < 2; i++)
	{
		MPoint scalePivot = i == 0 ? cvPoints[0] : cvPoints[cvPoints.size() - 1];
		bool shouldCreateCap = shouldCreateCaps[i];
		bool shouldFlipVert = shouldFlipVerts[i];
		int pipeVertexOffset = pipeVertexOffsets[i];
		int dataIndex = dataIndices[i];
		int extrusionDir = extrusionDirs[i];
		int face = faces[i];
		double twist = twists[i];
		MVector spanVectorDirection = (i == 0 ? spanVectors[0] : spanVectors[spanVectors.size() - 1]).normal();

		// back up existing verts a bit for edge loop
		int* vertexIndices = ArrayUtil::range(pipeVertexOffset, pipeVertexOffset + numSides);
		MTransformationMatrix transform;
		double capThick = shouldCreateCap ? capThickness : 0;
		MVector translateVector = spanVectorDirection * (-extrusionDir * (capThick + initialCapTightness));
		transform.setTranslation(translateVector, MSpace::kWorld);
		MeshUtil::transformPoints(fnMesh, vertexIndices, numSides, transform);
		delete[] vertexIndices;
		scalePivot += translateVector;

		
		// extrude down for initial edge loop 1
		MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, spanVectorDirection, extrusionDir * initialCapTightness);

		if (!shouldCreateCap) continue;

		// extrude out for initial edge loop 2
		MeshUtil::doExtrudeWithUniformScale(fnMesh, vertexOffset, numSides, face, scalePivot, (sectionRadius + initialCapTightness) / sectionRadius);
		
		// make cap shape; 0 = circle so do nothing
		if (numCapSides > 0)
		{
			int* vertexIndices = ArrayUtil::range(vertexOffset, vertexOffset + numSides);
			 
			// get around weird maya behavior
			if (shouldFlipVert)
			{
				std::swap(vertexIndices[0], vertexIndices[1]);
			}

			MTransformationMatrix capTransform;
			MVector vRotSrc(0, 0, 1);
			MVector vRotDst = spanVectorDirection * extrusionDir;
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
		MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, spanVectorDirection, extrusionDir * secondCapTightness);
		
		// body extrusion
		MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, spanVectorDirection, extrusionDir * (capThickness - secondCapTightness - thirdCapTightness));
		
		// extrude down for third edge loop 1
		MeshUtil::doExtrudeWithTranslate(fnMesh, vertexOffset, numSides, face, scalePivot, spanVectorDirection, extrusionDir * thirdCapTightness);

		// extrude in for third edge loop 2
		MeshUtil::doExtrudeWithUniformScale(fnMesh, vertexOffset, numSides, face, scalePivot, (capRadius - thirdCapTightness) / capRadius);
		
		// extrude in
		MeshUtil::doExtrudeWithUniformScale(fnMesh, vertexOffset, numSides, face, scalePivot, 0.5);

		// bolts TODO
	}
	return MS::kSuccess;
}

int PipeNode::numVertsInPipe()
{
	int numEdgeLoops = 2;
	for (int i = 0; i < angles.size(); i++)
	{
		numEdgeLoops += cornerDivisions[i] + 1;
	}
	return numEdgeLoops * numSides;
}

void* PipeNode::creator()
{
	return new PipeNode;
}

MStatus PipeNode::initialize()
{ 
	MStatus stat;
	MFnTypedAttribute typedAttr;
	MFnNumericAttribute numAttr;

	attrObjMesh = typedAttr.create("mesh", "m", MFnData::kMesh);
	typedAttr.setKeyable(false);
	typedAttr.setWritable(false);
	addAttribute(attrObjMesh);

	attrObjCurve = typedAttr.create("curve", "c", MFnData::kNurbsCurve);
	typedAttr.setKeyable(false);
	typedAttr.setReadable(false);
	addAttribute(attrObjCurve);
	attributeAffects(attrObjCurve, attrObjMesh);

	AttributeUtil::registerNumericAttribute<double>(&attrObjSectionRadius, &attrObjMesh, "radius", "r", MFnNumericData::kDouble, 1.0);
	AttributeUtil::registerNumericAttribute<int>(&attrObjRadialSubdivisions, &attrObjMesh, "radialSubdivisions", "d", MFnNumericData::kInt, 8, true, 3, true, 100);
	AttributeUtil::registerNumericAttribute<int>(&attrObjNumCapSides, &attrObjMesh, "capShape", "cs", MFnNumericData::kInt, 4);
	AttributeUtil::registerNumericAttribute<double>(&attrObjCapThickness, &attrObjMesh, "capThickness", "ct", MFnNumericData::kDouble, 1.0);
	AttributeUtil::registerNumericAttribute<bool>(&attrObjCreateBeginningCap, &attrObjMesh, "createBeginningCap", "cbc", MFnNumericData::kBoolean, true);
	AttributeUtil::registerNumericAttribute<double>(&attrObjBeginningCapRotation, &attrObjMesh, "beginningCapRotation", "bcr", MFnNumericData::kDouble, 0.0, true, 0.0, true, 360.0);
	AttributeUtil::registerNumericAttribute<bool>(&attrObjCreateEndCap, &attrObjMesh, "createEndCap", "cec", MFnNumericData::kBoolean, true);
	AttributeUtil::registerNumericAttribute<double>(&attrObjEndCapRotation, &attrObjMesh, "endCapRotation", "ecr", MFnNumericData::kDouble, 0.0, true, 0.0, true, 360.0);
	AttributeUtil::registerNumericAttribute<double>(&attrObjCapRadius, &attrObjMesh, "capRadius", "cra", MFnNumericData::kDouble, 2.0, true, 0.0);
	AttributeUtil::registerNumericAttribute<double>(&attrObjInitialCapTightness, &attrObjMesh, "initialCapTightness", "ict", MFnNumericData::kDouble, 0.1, true, 0.001);
	AttributeUtil::registerNumericAttribute<double>(&attrObjSecondCapTightness, &attrObjMesh, "secondCapTightness", "sct", MFnNumericData::kDouble, 0.1, true, 0.001);
	AttributeUtil::registerNumericAttribute<double>(&attrObjThirdCapTightness, &attrObjMesh, "thirdCapTightness", "tct", MFnNumericData::kDouble, 0.1, true, 0.001);
	AttributeUtil::registerNumericAttribute<double>(&attrObjSideCapTightness, &attrObjMesh, "sideCapTightness", "sict", MFnNumericData::kDouble, 0.1, true, 0.001);
	
	attrObjRadiusArray = numAttr.create("radiusPerCorner", "rpc", MFnNumericData::kDouble, 2.0);
	numAttr.setArray(true);
	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	addAttribute(attrObjRadiusArray);
	attributeAffects(attrObjRadiusArray, attrObjMesh);

	AttributeUtil::registerNumericAttribute<double>(&attrObjPipeRotation, &attrObjMesh, "pipeRotation", "pr", MFnNumericData::kDouble, 0.0, true, 0.0, true, 360.0);
	AttributeUtil::registerNumericAttribute<double>(&attrObjCornerAnglePerDivision, &attrObjMesh, "cornerAnglePerDivision", "cad", MFnNumericData::kDouble, 30.0, true, 0.1, true, 90.0);
	
	attrObjBoltPath = typedAttr.create("boltPath", "bp", MFnData::kString, MFnStringData().create());
	typedAttr.setKeyable(true);
	typedAttr.setStorable(true);
	addAttribute(attrObjBoltPath);
	attributeAffects(attrObjBoltPath, attrObjMesh);

	AttributeUtil::registerNumericAttribute<double>(&attrObjBoltScale, &attrObjMesh, "boltScale", "bs", MFnNumericData::kDouble, 1.0);
	AttributeUtil::registerNumericAttribute<int>(&attrObjSeed, &attrObjMesh, "seed", "s", MFnNumericData::kInt, 2938);
	return MStatus::kSuccess;
}
