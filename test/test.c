#include <stdio.h>
#include <malloc.h>

int testarray(double in[],int n, double out[])
{
    int i;
    for (i = 0; i < n; i++)
    {
        out[i] = in[i]*2.0;// clone
        int j=1;
    }
    return 1;
}

int main()
{
    double kk[1000000];
    //kk = (double*)malloc(1000000*sizeof(double));
    for (int j = 0; j < 10000000; j++)
    {
        kk[j] = j;
    }
    int sz = (int)(sizeof(kk)/sizeof(kk[0]));
    double* out = (double*)malloc(sz*sizeof(double));
        testarray(kk,sz,out);
    for (int i = 0; i < sz;i++)
    {
        printf("%lf    \n", *(out+i));
    }
}