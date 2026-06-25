// Array Manipulation by Field Name

#include "parameters.h"
#include "string.h"
#include "mex.h"


// 
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE
#define error(name) mexPrintf("Not Found :%s\n",name)
#endif
#else
	#include <assert.h>
	#define error(name) printf("Not Found :%s\n",name) ; assert(false)
#endif

//#include <assert.h>

// Find index corresponding to name (parameter name) in the field paramfield

int Find_Index(const char* name, parameters * paramfield, int tab_size)
{
	int ret = -1;
    int i ;
	bool found = false;
	for (i = 0; i < tab_size && found ==false; i++)
	{
		if (strcmp(paramfield[i].ParamName, name) == 0)
		{
			ret = i;
			found = true;
				
		}
	}
	if (!found)
	{
		//error(name) ;
		ret = -1;
	}
	return ret;
}
int getIntValue(double * tab, const char * field_name, parameters * paramfield, int tab_size)
{
	int index = Find_Index(field_name,paramfield,tab_size);
    if (index >=0)
	{
		return (int)tab[index];
	}
	return 0;

}

bool getBoolValue(double * tab, const char * field_name, parameters * paramfield, int tab_size)
{
	int index = Find_Index(field_name,paramfield,tab_size);
    if (index >=0)
	{
		return (bool)tab[index];
	}
	return false;

}

float getFloatValue(double * tab, const char * field_name, parameters * paramfield, int tab_size)
{
	int index = Find_Index(field_name,paramfield,tab_size);
    if (index >=0)
	{
		return (float)tab[index];
	}
	return 0.0f;

}

double getDoubleValue(double * tab, const char * field_name, parameters * paramfield, int tab_size)
{
	int index = Find_Index(field_name,paramfield,tab_size);
    if (index >=0)
	{
		return tab[index];
	}
	return 0.0;
}

void SetValue(double * tab, const char * field_name, parameters * paramfield, int tab_size, double value)
{
	int index = Find_Index(field_name,paramfield,tab_size);
    if (index >=0)
	{
		tab[index] = value;
	}
}