#include <stdio.h>
#include "msis00.h"
int main()
{
	int day = 172;
	float sec = 29000.;
	float ALT=50.;
	float lon=-70;
	float lat=60;
	float t = gettemp(day,sec,ALT,lat,lon);
	printf("%.2f\n",t);
	return 0;
	
}
