#include <math.h>
const  double   dia= 6378137.0    ;                      
const  double   fai1= 1/298.257223563  ;                 
const  double   pai= 3.14159265358979 ; 
const  double 	e12=0.00669438499959;    



void  DBLH_DXYZ(double B,double L,double H,double *aX,double *aY,double *aZ)
{
	double di_N;
	double W;
	double X,Y,Z;
	W=sqrt(1.0-e12*sin(B)*sin(B));
	di_N=dia/W;
	X=(di_N+H)*cos(B)*cos(L);
	Y=(di_N+H)*cos(B)*sin(L);
	Z=(di_N*(1-e12)+H)*sin(B); 
	*aX = X;
	*aY = Y;
	*aZ = Z;
}

//
void DXYZ_DBLH(double X,double Y,double Z,double *aB,double *aL,double *aH)
{
	double l;
	double a;
	double temp;
	double s1,s2;
	double B,L,H;
	//�����ؾ���
	l=atan(Y/X);
	if((X>0)&&(Y>0))
	{
		L=l;   
	}
	else if((X<0)&&(Y>0))
	{
		L=l+pai;
	}
	else if((X<0)&&(Y<0))
	{
		L=l-pai;
	}
	else 
	{
		L=l;
	}
	
	a=atan(Z*(1.0+e12)/sqrt(X*X+Y*Y));
	
	do
	{
		temp=a;
		s1=dia*e12*tan(a) ;
		s2=sqrt(1.0+(1.0-e12)*tan(a)*tan(a));
		a=atan(( Z+(s1/s2))/sqrt(X*X+Y*Y));	
		
	}while((a-temp)>=0.0000000000001);
	
    B=a;
	
	H=Z/sin(B)-(dia/sqrt(1-e12*sin(B)*sin(B)))*(1.0-e12);
	*aB = B;
	*aL = L;
	*aH = H;
}