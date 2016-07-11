#include <iostream>
#include <cstdio>
#include <cmath>
#include "FdtdManager.h"

int main()
{
	FdtdManager manager;
	manager.init("BoreHole.json");
	//manager.loadJson((char*)"a.json");
	//manager.buildModel();
	manager.beginFdtd();
	return 0;
}