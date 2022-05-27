#include "pch.h"

#include "TestVector.h"

bool assertEquals(std::vector<int> v1, std::vector<int> v2)
{
	bool assertEqual = true;

	for (int i = 0; i < v1.size(); i++)
	{
		assertEqual = assertEqual && (v1[i] == v2[i]);
	}

	return assertEqual;
}

bool assertEquals(std::vector<short> v1, std::vector<short> v2)
{
	bool assertEqual = true;

	for (int i = 0; i < v1.size(); i++)
	{
		assertEqual = assertEqual && (v1[i] == v2[i]);
	}

	return assertEqual;
}

bool assertEquals(std::vector<double> v1, std::vector<double> v2)
{
	bool assertEqual = true;

	for (int i = 0; i < v1.size(); i++)
	{
		assertEqual = assertEqual && (v1[i] == v2[i]);
	}

	return assertEqual;
}