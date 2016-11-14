#pragma once

#include <maya/MPxCommand.h>
#include <maya/MObject.h>
#include <maya/MDagModifier.h>
#include <maya/MSelectionList.h>

#include <vector>

class PipeGeneratorCommand : public MPxCommand
{
public:
	PipeGeneratorCommand();
	~PipeGeneratorCommand();
	MStatus doIt(const MArgList& args);
	MStatus redoIt();
	MStatus undoIt();
	static void* creator();
	bool isUndoable() const;
private:
	std::vector<MObject> curveNodeObjs;
	std::vector<MDagModifier*> dagModifiers;
	std::vector<MDGModifier*> dgModifiers;
	std::vector<MObject> transformNodeObjs;
	std::vector<MObject> shapeNodeObjs;
	std::vector<MObject> pipeNodeObjs;
	std::vector<int> instanceNumbers;
	MSelectionList selectionList;

	static int instanceCounter;
};