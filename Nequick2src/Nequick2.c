#include "Nequick2.h"
extern double nequick_(double *h,double *alat, double *along, int *ath,double *flx, double *UT);
double runNequick2(double h,double alat,double along,int month,double flx,double UT)
{
	return nequick_(&h,&alat,&along,&month,&flx,&UT);
}
