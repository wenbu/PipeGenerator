#pragma once

#include <maya/MPoint.h>
#include <maya/MTransformationMatrix.h>

namespace TransformUtil
{
	MTransformationMatrix getUniformScale(double s, const MPoint& pivot);
}