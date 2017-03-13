#pragma once

#include <maya/MObject.h>
#include <maya/MFnNumericData.h>
#include <maya/MStatus.h>

namespace AttributeUtil
{
	// probably a better way to do this than pass in bools
	// assumes exactly one affected object
	template<typename T> static MStatus registerNumericAttribute(MObject* attribute, MObject* affectedObject, const char* fullName, const char* shortName, MFnNumericData::Type type,
		T defaultValue = 0, bool hasMin = false, T minValue = 0, bool hasMax = false, T maxValue = 0);

	template<typename T>
	MStatus registerNumericAttribute(MObject* attribute, MObject* affectedObject, const char * fullName, const char * shortName, MFnNumericData::Type type, T defaultValue, bool hasMin, T minValue, bool hasMax, T maxValue)
	{
		MStatus stat;
		MFnNumericAttribute numAttr;
		*attribute = numAttr.create(fullName, shortName, type, defaultValue, &stat);
		if (!stat) return stat;
		stat = numAttr.setKeyable(true);
		if (!stat) return stat;
		stat = numAttr.setReadable(true);
		if (!stat) return stat;
		if (hasMin)
		{
			stat = numAttr.setMin(minValue);
			if (!stat) return stat;
		}
		if (hasMax)
		{
			stat = numAttr.setMax(maxValue);
			if (!stat) return stat;
		}
		stat = MPxNode::addAttribute(*attribute);
		if (!stat) return stat;
		stat = MPxNode::attributeAffects(*attribute, *affectedObject);
		if (!stat) return stat;
		return MS::kSuccess;
	}
}