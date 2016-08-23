#include "fintrf.h"
#include "NeQuick_2.for"
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
         IMPLICIT REAL*8 (A-H, O-Z)
         integer nlhs, nrhs
         mwPointer plhs(*), prhs(*)
         mwPointer mxGetPr, mxGetM, mxGetN
         mwPointer mxCreateDoubleMatrix
         real*8 NeQuick
         real*8 h,alat,along,flx,UT,y,mth
         mwPointer m,n,x_pr,y_pr
         mwSize size
         x_pr = mxGetPr(prhs(1))
         call mxCopyPtrToReal8(x_pr,h,1)
C         h = data(1)
C         alat=data(2)
C         along=data(3)
C         mth=data(4)
C         flx=data(5)
C         UT=data(6)
         x_pr = mxGetPr(prhs(2))
         call mxCopyPtrToReal8(x_pr,alat,1)
         x_pr = mxGetPr(prhs(3))
         call mxCopyPtrToReal8(x_pr,along,1)
         x_pr = mxGetPr(prhs(4))
         call mxCopyPtrToReal8(x_pr,mth,1)
         x_pr = mxGetPr(prhs(5))
         call mxCopyPtrToReal8(x_pr,flx,1)
         x_pr = mxGetPr(prhs(6))
         call mxCopyPtrToReal8(x_pr,UT,1)
         plhs(1) = mxCreateDoubleMatrix(1,1,0)
C         plhs(1) = mxCreateNumericMatrix(1, 1,
C     +                mxClassIDFromClassName('int16'), 0)
C         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput',
C     +                           'One input required.')
         y_pr = mxGetPr(plhs(1))
C         y = mth
         y = NeQuick(h,alat,along,mth,flx,UT)

         call mxCopyReal8ToPtr(y,y_pr,1)

C         call mexErrMsgIdAndTxt('MATLAB','Debug')
      end subroutine
