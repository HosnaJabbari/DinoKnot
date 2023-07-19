#ifndef HFOLD_INTERACTING_H_
#define HFOLD_INTERACTING_H_

#include "Result.h"
#include <vector>


int validateModelType(char* type);  
double get_START_HYBRID_PENALTY(int type1, int type2);
int isDirectory(const char *path);

#endif /*HFOLD_INTERACTING_H_*/
