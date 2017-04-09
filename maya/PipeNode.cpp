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
MObject PipeNode::attrObjPipeRotation;
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

	MDataHandle dataHandlePipeRotation = dataBlock.inputValue(attrObjPipeRotation, &stat); if (!stat) return stat;

	MDataHandle dataHandleOutMesh = dataBlock.outputValue(attrObjMesh, &stat); if (!stat) return stat;

	MObject curve = dataHandleInCurve.asNurbsCurveTransformed();
	numSides = dataHandleInRadialSubdivisions.asInt();
	sectionRadius = dataHandleSectionRadius.asDouble();
	cornerAnglePerDivision = dataHandleCornerAnglePerDivision.asDouble() * (2 * M_PI / 360);
	pipeRotation = dataHandlePipeRotation.asDouble() * (2 * M_PI / 360);
	
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

	outputData = dataCreator.create(&stat);
	if (!stat) return stat;

	stat = getPipeData();
	if (!stat) return stat;
	stat = makeBaseMesh();
	if (!stat) return stat;
	stat = makePipe();
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

	attrObjRadiusArray = numAttr.create("radiusPerCorner", "rpc", MFnNumericData::kDouble, 2.0);
	numAttr.setArray(true);
	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	addAttribute(attrObjRadiusArray);
	attributeAffects(attrObjRadiusArray, attrObjMesh);

	AttributeUtil::registerNumericAttribute<double>(&attrObjPipeRotation, &attrObjMesh, "pipeRotation", "pr", MFnNumericData::kDouble, 0.0, true, 0.0, true, 360.0);
	AttributeUtil::registerNumericAttribute<double>(&attrObjCornerAnglePerDivision, &attrObjMesh, "cornerAnglePerDivision", "cad", MFnNumericData::kDouble, 30.0, true, 0.1, true, 90.0);
	
	return MStatus::kSuccess;
}
