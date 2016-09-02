#include <stdio.h>
#include "xyz2blh.h"
#define PI 3.14159265
int main()
{
    double b = 35.4/180*PI;
    double l = -110.0/180*PI;
    double h = 30*1000;
    double x,y,z;
    double B,L,H;
    DBLH_DXYZ(b,l,h,&x,&y,&z);
    printf("%lf\t%lf\t%lf\n",x,y,z);
    DXYZ_DBLH(x,y,z,&B,&L,&H);
    printf("%lf\t%lf\t%lf\n",B/PI*180,L/PI*180,H);
    return 0;
}