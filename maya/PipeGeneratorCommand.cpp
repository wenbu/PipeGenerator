#include "PipeGeneratorCommand.h"
#include "PipeNode.h"

#include <maya/MGlobal.h>
#include <maya/MItSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnTransform.h>

PipeGeneratorCommand::PipeGeneratorCommand()
{
}


PipeGeneratorCommand::~PipeGeneratorCommand()
{
	// need to delete MDGModifiers and MDagModifiers?
}

int PipeGeneratorCommand::instanceCounter = 0;

MStatus PipeGeneratorCommand::doIt(const MArgList& args)
{
	MStatus stat;
	selectionList = MSelectionList();
	stat = MGlobal::getActiveSelectionList(selectionList);
	if (!stat)
	{
		stat.perror("Failed to get active selection list");
		return stat;
	}
	
	return redoIt();
}

MStatus PipeGeneratorCommand::redoIt()
{
	MStatus stat;
	MItSelectionList selectionListIterator = MItSelectionList(selectionList);

	while (!selectionListIterator.isDone())
	{
		// get dagpath to selected item
		MDagPath dagPath;
		stat = selectionListIterator.getDagPath(dagPath);
		if (!stat)
		{
			stat.perror("Failed to get dagpath from selection list iterator");
			return stat;
		}

		// get shape
		stat = dagPath.extendToShape();
		if (!stat)
		{
			stat.perror("Failed to extend dagpath to shape");
			return stat;
		}

		MObject nodeObj = dagPath.node(&stat);
		if (!stat)
		{
			stat.perror("Failed to get node MObject from dagpath");
			return stat;
		}

		if (nodeObj.hasFn(MFn::kNurbsCurve))
		{
			curveNodeObjs.push_back(nodeObj);
			MDagModifier* dagModifier = new MDagModifier;
			MDGModifier* dgModifier = new MDGModifier;

			// create nodes
			MObject meshTransformObj = dagModifier->createNode("transform", MObject::kNullObj, &stat);
			if (!stat)
			{
				stat.perror("Failed to create transform node");
				return stat;
			}
			MObject meshShapeObj = dagModifier->createNode("mesh", meshTransformObj, &stat);
			if (!stat)
			{
				stat.perror("Failed to create shape node");
				return stat;
			}
			MObject pipeNodeObj = dgModifier->createNode(PipeNode::nodeID, &stat);
			if (!stat)
			{
				stat.perror("Failed to create pipe node");
				return stat;
			}

			dagModifiers.push_back(dagModifier);
			dgModifiers.push_back(dgModifier);
			transformNodeObjs.push_back(meshTransformObj);
			shapeNodeObjs.push_back(meshShapeObj);
			pipeNodeObjs.push_back(pipeNodeObj);
			instanceNumbers.push_back(PipeGeneratorCommand::instanceCounter);
			instanceCounter++;
		}

		stat = selectionListIterator.next();
		if (!stat)
		{
			stat.perror("Failed to get next item from selection list iterator");
			return stat;
		}
	}

	// assume all vectors are same size
	// if not then something is wrong lol!!
	for (int i = 0; i < curveNodeObjs.size(); i++)
	{
		stat = dagModifiers[i]->doIt();
		if (!stat)
		{
			stat.perror("dagModifier failed to do it");
			return stat;
		}
		stat = dgModifiers[i]->doIt();
		if (!stat)
		{
			stat.perror("dgModifier failed to do it");
			return stat;
		}

		MFnDependencyNode pipeDepNodeFn(pipeNodeObjs[i]);
		MFnDependencyNode curveDepNodeFn(curveNodeObjs[i]);
		MFnDependencyNode shapeDepNodeFn(shapeNodeObjs[i]);

		// get plugs
		MPlug curveSrcPlug = curveDepNodeFn.findPlug("worldSpace", &stat);
		if (!stat)
		{
			stat.perror("Failed to get worldSpace plug from curve");
			return stat;
		}
		curveSrcPlug = curveSrcPlug.elementByLogicalIndex(0, &stat);
		if (!stat)
		{
			stat.perror("Failed to get worldSpace[0] plug");
			return stat;
		}
		MPlug curveDstPlug = pipeDepNodeFn.findPlug("curve", &stat);
		if (!stat)
		{
			stat.perror("Failed to get curve plug from pipe node");
			return stat;
		}
		MPlug meshSrcPlug = pipeDepNodeFn.findPlug("mesh", &stat);
		if (!stat)
		{
			stat.perror("Failed to get mesh plug from pipe node");
			return stat;
		}
		MPlug meshDstPlug = shapeDepNodeFn.findPlug("inMesh", &stat);
		if (!stat)
		{
			stat.perror("Failed to get inMesh plug from shape node");
			return stat;
		}

		MFnTransform transformFn(transformNodeObjs[i]);
		transformFn.setName("coolPipe" + instanceNumbers[i]);
		pipeDepNodeFn.setName("coolPipeMaker" + instanceNumbers[i]);
		shapeDepNodeFn.setName("coolPipeShape" + instanceNumbers[i]);

		// connect
		MDGModifier connDgModifier;
		stat = connDgModifier.connect(curveSrcPlug, curveDstPlug);
		if (!stat)
		{
			stat.perror("Failed to connect curve plugs");
			return stat;
		}
		stat = connDgModifier.connect(meshSrcPlug, meshDstPlug);
		if (!stat)
		{
			stat.perror("Failed to connect mesh plugs");
			return stat;
		}
		stat = connDgModifier.commandToExecute("select " + transformFn.name());
		if (!stat)
		{
			stat.perror("Failed to execute command to select transform node");
			return stat;
		}
		stat = connDgModifier.commandToExecute("hyperShade -assign lambert1");
		if (!stat)
		{
			stat.perror("Failed to execute command to assign lambert1");
			return stat;
		}
		stat = connDgModifier.doIt();
		if (!stat)
		{
			stat.perror("connDgModifier failed to do it");
			return stat;
		}

		// copy pipe mesh data to outMesh on shape node, since update isn't immediate
		MPlug shapeOutMeshPlug = shapeDepNodeFn.findPlug("outMesh", &stat);
		if (!stat)
		{
			stat.perror("Failed to find outMesh plug on shape node");
			return stat;
		}
		MObject meshDataObj = meshSrcPlug.asMObject(MDGContext::fsNormal, &stat);
		if (!stat)
		{
			stat.perror("Failed to get mesh data from pipe node");
			return stat;
		}
		stat = shapeOutMeshPlug.setMObject(meshDataObj);
		if (!stat)
		{
			stat.perror("Failed to set mesh data on shape node");
			return stat;
		}
		// need try-catch?
	}
	return MS::kSuccess;
}

MStatus PipeGeneratorCommand::undoIt()
{
	for (int i = 0; i < dagModifiers.size(); i++)
	{
		dgModifiers[i]->undoIt();
		dagModifiers[i]->undoIt();
	}
	return MS::kSuccess;
}

void* PipeGeneratorCommand::creator()
{
	return new PipeGeneratorCommand;
}

bool PipeGeneratorCommand::isUndoable() const
{
	return true;
}