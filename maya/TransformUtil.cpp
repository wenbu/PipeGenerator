#include "TransformUtil.h"

MTransformationMatrix TransformUtil::getUniformScale(double s, const MPoint& pivot)
{
	MTransformationMatrix transform;
	double3 scale = { s, s, s };
	transform.setScale(scale, MSpace::kWorld);
	transform.setScalePivot(pivot, MSpace::kWorld, false);
	return transform;
}
