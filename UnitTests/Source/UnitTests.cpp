// UnitTests.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "UnitTests.h"

#include "TestBinaryDecimalConversion.h"

// TODO: This is an example of a library function
bool RunUnitTests()
{
	std::vector<double> v1 = { 0.0, 0.0, 0.0 };
	std::vector<double> v2 = { 1.0, 0.0, 0.0 };

	if (assertEquals(v1, v2))
		std::cout << "I called the library function!\n";

	TestMinConversions();

	return true;
}
