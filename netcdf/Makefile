testcdf: testnetcdfmain.o
	cc -o testcdf testnetcdfmain.o -L/usr/local/netcdf/lib -static -lnetcdf
testnetcdfmain.o: testnetcdfmain.c
	cc -c -I/usr/local/netcdf/include testnetcdfmain.c
