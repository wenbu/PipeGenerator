#include <iostream>
#include <maya/MSimple.h>
#include <maya/MIOStream.h>

DeclareSimpleCommand(HelloWorld, "Autodesk", "2016");

MStatus HelloWorld::doIt(const MArgList&)
{
	std::cout << "yo what's up\n" << std::endl;
	return MS::kSuccess;
}