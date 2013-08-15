
#ifndef _SYNC_H_
#define _SYNC_H_

#include "rate.h"
#include "cvector.h"
#include <stdio.h>
#include <math.h>
#include "utilities.h"
#include "const.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;

void LoadSyncTable(dCVector* syncTable, string aDirTables);
void InitializeSynchrotron(const double B_0, const dCVector* pEnergy, 
			   const dCVector* pEnergyWidth,
			   dCVector* synchrotronLoss, DiffRate* syncRate,
			   string aDirTables);

#endif
