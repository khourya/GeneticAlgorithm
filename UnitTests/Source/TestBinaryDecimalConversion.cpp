#include "pch.h"

#include "TestVector.h"
#include "TestBinaryDecimalConversion.h"

bool TestMinConversions()
{
	bool testsPass = true;

	testsPass = testsPass && TestMinConversions_1Variable();
	testsPass = testsPass && TestMinConversions_2SameLengthVariables();
	testsPass = testsPass && TestMinConversions_2MixedLengthVariables();

	return testsPass;
}

bool TestMinConversions_1Variable()
{
	int variableLength = 3;
	double minVal = 0.0;
	double maxVal = 10.0;

	DesignVariable* v1 = new DesignVariable(variableLength, minVal, maxVal);
	std::vector<DesignVariable*> designVars = { v1 };
	std::vector<double> vars = { minVal };

	Individual* indy = new Individual(designVars, vars);
	std::vector<short> actualChromosome = indy->GetChromosomeVector();

	std::vector<short> expectedChromosome(variableLength, 0);

	return assertEquals(actualChromosome, expectedChromosome);
}

bool TestMinConversions_2SameLengthVariables()
{
	return true;
}

bool TestMinConversions_2MixedLengthVariables()
{
	return true;
}

bool TestMaxConversions()
{
	return true;
}