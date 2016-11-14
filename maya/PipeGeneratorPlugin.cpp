#include "PipeNode.h"
#include "PipeGeneratorCommand.h"

#include <maya/MObject.h>
#include <maya/MFnPlugin.h>

MStatus initializePlugin(MObject obj)
{
	MFnPlugin plugin(obj, "bwu", "1.0", "Any");
	MStatus stat = plugin.registerNode(PipeNode::nodeName, PipeNode::nodeID, PipeNode::creator, PipeNode::initialize, MPxNode::kDependNode, &PipeNode::nodeClassify);
	if (!stat)
	{
		stat.perror("Failed to register PipeNode");
		return stat;
	}
	stat = plugin.registerCommand("createPipeNode", PipeGeneratorCommand::creator);
	if (!stat)
	{
		stat.perror("Failed to register createPipeNode command");
		return stat;
	}
	return MS::kSuccess;
}

MStatus uninitializePlugin(MObject obj)
{
	MFnPlugin plugin(obj);
	MStatus stat = plugin.deregisterNode(PipeNode::nodeID);
	if (!stat)
	{
		stat.perror("Failed to deregister PipeNode");
		return stat;
	}
	stat = plugin.deregisterCommand("createPipeNode");
	if (!stat)
	{
		stat.perror("Failed to deregister createPipeNode command");
		return stat;
	}
	return MS::kSuccess;
}