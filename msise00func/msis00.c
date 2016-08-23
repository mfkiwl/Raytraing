#include "msis00.h"
#include <stdio.h>
extern void gtd7_(int *IYD, float*SEC, float*ALT, float *GLAT, float *GLONG, float *STL, float *F107A, float *F107, float *AP, int *MASS, float *D, float *T);
extern void GHP7_(int *IYD, double *SEC, double *ALT, double *GLAT, double *GLONG, double *STL, double *F107A, double *F107, double *AP, double *D, double *T, double PRES);
float gettemp(int iday, float sec, float altkm, float glatdeg, float glondeg)
{
	float stl = sec/3600. + glondeg/15.;
	stl = 16;
	float f107=150.;
	float f107a=150.;
	float ap[7] = {4.};
	int mass = 0;
	float d[9]={0};
	float t[2]={0};
	gtd7_(&iday,&sec,&altkm,&glatdeg,&glondeg,&stl,&f107a,&f107,ap,&mass,d,t);
	return t[1];
}
