#include "stdafx.h"
#include "Grid.h"
#include <windows.h>
#include <iostream>
#include <fstream>

int main()
{
	Grid test;
	test.generateGrid("");  //enter data file's path
	test.aggregateGrid();
	test.calculateTemperatures();
	Sleep(5000);
}