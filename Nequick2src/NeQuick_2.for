C =========================================================================
C                  NeQuick2 P.531 electron density model
C     		Release date 22 Oct 2011
C     		(Internal developers reference 2.1.0)
C     ========================================================================
C     DISCLAIMER:
C	This software is meant for scientific use only.
C         Please acknowledge the Aeronomy and Radiopropagation Laboratory
C     	of the Abdus Salam International Centre for Theoretical Physics
C     	Trieste, Italy
C
C     This software is provided to The User "as is". ITU-R and the authors assume 
C	no responsibility whatsoever for its use by other parties, and makes no 
C	guarantees, expressed or implied, about the quality, reliability, or any 
C 	other characteristic of this software. Under no circumstances shall ITU, 
C	the author(s) or contributor(s) be liable for damages resulting directly or 
C	indirectly from the use, misuse or inability to use this software.
C	
C    The user may not modify this software in content and then present it or 
C	results derived from it as ITU-R or author(s) material. 
C     =========================================================================
C
C  NeQuick was developed at 
C  the Abdus Salam International Centre for Theoretical
C  Physics, Trieste, Italy
C  and at the University of Graz, Austria.
C
C NeQuick model references:
C Radicella, S.M., R. Leitinger, The evolution of th DGR approach to
C model electron density profiles, Adv. Space Res., 27(1), 35-40, 2001.
C
C Hochegger, G., B. Nava, S.M. Radicella, R. Leitinger, A family of
C ionospheric models for different uses- Phys. Chem. Earth (C), 
C 25 (4), 307-310, 2000.
C
C Leitinger, R., M.L. Zhang, S.M. Radicella, An improved bottomside
C for the ionospheric electron density model NeQuick, Annals of
C Geophysics 48(3) 525-534, 2005.
C
C Coisson, P., S.M. Radicella, R. Leitinger, B. Nava, Topside electron 
C density in IRI and NeQuick: features and limitation, 
C Adv. Space Res. 37, 937-942, 2006.
C
C NeQuick is based on the Di Giovanni - Radicella (DGR) model
C   (Di Giovanni, G., S.M. Radicella, An analytical model of the electron
C      density profile in the ionosphere, Adv. Space Res. 10, 27-30, 1990)
C   which was modified to the requirements of the COST 238 Action PRIME
C   to give vertical electron content from ground to 1000 km consistent
C   with the COST 238 regional electron content model
C   (Radicella, S.M., M.-L. Zhang, The improved DGR analytical model of
C      electron density height profile and total electron content in the
C      ionosphere. A. Geofisica, 38, 35-41, 1995).
C
C The topside of NeQuick is a simplified approximation to a diffusive
C   equilibrium, the main improvement over the DGR and COST 238 models
C   being a limited increase with height of the electron density scale
C   height used.
C
C  NeQuick is a "profiler" which makes use of three profile anchor points:
C   E layer peak (at a fixed height of 120 km), F1 peak, F2 peak.
C   To model the anchor points it uses the "ionosonde parameters"
C   foE, foF1, foF2 (critical frequencies) and M(3000)F2 (transfer
C   parameter). For foE it uses a model by John Titheridge; foF1 is taken
C   to be proportional to foE during daytime (foF1=1.4*foE) and 0 during
C   nighttime. For foF2 and M(3000)F2 it uses the ITU-R (CCIR)
C   coefficients in the same way used by the International Reference
C   Ionosphere (IRI).
C
C The bottom side of the electron density profile consists of the
C   superposition of three Epstein layers which peak at the anchor points.
C   The Epstein layers have different thickness parameters for their
C   bottom and top sides (5 "semi-Epstein" layers). The topside of the
C   electron density profile consists of the topside of an Epstein
C   layer with a height dependent thickness parameter.
C
C The sub-models contained in NeQuick use monthly average values of solar
C   activity in two forms: average sunspot number R12 and 10.7 cm
C   solar radio flux F107. The latter is considered to be the primary
C   input parameter. A fixed relation between R12 and F107 is used:
C            F107= 63.7+8.9D-4*R12**2+0.728*R12  (ITU-R recommendation)
C       or   R12 = sqrt(167273.0+(flx-63.7)*1123.6)-408.99
C   NeQuick calculates modip by third order interpolation in
C   geographic latitude and longitude of the grid contained in the file
C   modip.asc.
C
C Authors of the modified Di Giovanni - Radicella (DGR) model code:
C   Man-Lian Zhang and S.M. Radicella, Abdus Salam ICTP Trieste.
C
C   Source of the code for the ITU-R coefficients (subroutines cciri and
C   function gamma1) and for the   ccirXX.asc   data files:
C   International Reference Ionosphere (IRI)
C   (  real FUNCTION xNeIRI(height,lati,longi,dhour,month,iday,R)
C      SUBROUTINE F2OUT(XMODIP,XLATI,XLONGI,FF0,XM0,UT,
C     &                              FOF2,XM3000)
C      REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,IHARM,NQ,
C    in:
C    IRIF12.FOR, Version 12.2.2, September 1994  )
C
C   Author of modifications of DGR code,
C   of adaptations of CCIR coefficients code
C   and author of the original NeQuick package:
C    Dr. Reinhart Leitinger.
C   the Abdus Salam International Centre for Theoretical Physics
C   strada costiera 11
C   34014 Trieste (TS)
C   Italy
C   e-mail: bnava@ictp.it
C           yenca@ictp.it
C           
C
C   with contributions by Bruno Nava,Gerald Hochegger,Johanna Wagner,
C   Pierdavide Coisson and Yenca Migoya Orue'
C
C =========================================================================
C
C This file is written in FORTRAN 77.
C      Implicit double precision for all real variables and
C      functions (real*8). It consists of 410 FORTRAN statement lines.
C Driver programs delivered:
C      source programs:
C	slQu_2.for    usage: slQu_2
C			provides a text-based interface step-by-step procedure to calculate
C			electron density and total electron content along any 				            
C			arbitrarily chosen straight line from a number of input parameters.
C			The output file slQu.dat is created, which contains the electron 
C			density in units of [m-3] along the profile (if required) and TEC 
C			in units 1015 [m-2].
C	NeQVal.for 	usage: NeQVal F10.7 InputFile (> OutputFile) 
C			calculates Slant Total Electron Content for a given Solar Flux at 
C			10.7 cm (F10.7) all lines in InputFile, where 
C			each line defines the epoch and the coordinates of begin and end 
C			points, the year. The results is provided in the Standard Output 
C			with a line including the same line of the input file plus an addtional 
C			column with the Slant Total Electron Content along the ray.
C
C =========================================================================
C
C program modules
C      real*8 function NeQuick(h,alat,along,mth,flx,UT)
C      subroutine prepNeQ(alat,along,mth,UT,flx,hm,aep,bb)
C      function*8 NeNeQ(h,hm,aep,bb)
C      real*8 function NeMdGR(aep,hm,bb,h)
C      real*8 function topq(h,No,hmax,Ho)
C      subroutine prepmdgr(mth,R12,foF2,efoF1,efoE,M3000,hm,bb,aep)
C      subroutine ef1(alat,mth,flx,chi,foF2,efoE,efoF1)
C      subroutine cciri(xMODIP,mth,UT,R12,alat,along,foF2,M3000)
C      real*8 function gamma1(xMODIP,alat,along,hour,iharm,nq,
C     ,                          k1,m,mm,m3,sfe)
C      real*8 function peakh(foE,foF2,M3000)
C      subroutine sdec(mth,UT,sdelta,cdelta)
C      subroutine modin(filenam,pmodip)
C      real*8 function amodip(pmodip,alat,along)
C      real*8 function fexp(a)
C      real*8 function djoin(f1,f2,alpha,x)
C      real*8 function finter3(z,x)
C
C =========================================================================
C
C Data files needed:
C   (1) 12 CCIR coefficients files       ccir11.asc   (for January) to
C                                        ccir22.asc   (for December)
C   (2) modip grid file  modip.asc
C
C =========================================================================
C
C Input parameters:
C   NeQuick:
C     h:     height (km)
C     alat:  gg. latitude  (degrees N)
C     along: gg. longitude (degrees E)
C     mth:   month (1 .. 12)
C     flx:   10.7 cm solar radio flux (flux units)
C     UT:    Universal Time (hours)
C
C =========================================================================
C
C Output (function value) is electron density in units of m^-3
C   (electrons per cubic meter)
C
C =========================================================================
C
C
C Changes:
C 12.11.2007 optimization of prepmdgr and NeMdGR
C 14.11.2007 all constants in double precision
C 05.02.2010 control on a2 parameter to be non-negative
C 15.06.2010 initial value -1 instead of 0 in prepNeQ for flx0
C 21.06.2010 introduced modip file name as argument in subroutine modin


      real*8 function NeQuick(h,alat,along,mth,flx,UT)
      implicit real*8 (a-h,o-z)
      real*8 NeNeQ
      dimension hm(3),aep(3),bb(6)
      call prepNeQ(alat,along,mth,UT,flx,hm,aep,bb)
      NeQuick= NeNeQ(h,hm,aep,bb)
      return
      end

      subroutine prepNeQ(alat,along,mth,UT,flx,hm,aep,bb)
      implicit real*8 (a-h,o-z)
      real*8 M3000
      character*15 filenam
      dimension pmodip(0:183,-1:182)
      dimension hm(3),aep(3),bb(6)

      save

      parameter (pi12=0.261799387788149D0)
      parameter (DR=1.74532925199433D-2,RD=5.729577951308232D1)
      data UT0,mth0,flx0,jdip/-100.0D0,-1,-1.0D0,0/
      mth1=mth
      flx1=flx
      ut1=dmod(UT+24.0D0,24.0D0)

c     *** flux saturation above 193 FU and blocked below 63 FU to avoid
c     unrealistic or undefined electron density values! ***
      if (flx1 .gt. 193.0D0) then
      flx1=193.0D0
      write(*,'(2A/A/A)')'***WARNING! Solar flux limit F=193 (R12=150)',
     & 		' exceeded.',
     & 		' NeQuick saturates F to 193 units',
     & 		' (ITU-R P.1239 reccomendation).'
      endif
      if (flx1 .lt. 63.0D0) then
      write(*,'(2A/A)')'***WARNING! Solar flux below 63 FU (R12 <0)',
     &                   'program stopped!'
      stop
      endif

      if (jdip.eq.0) then
         filenam='modip.asc'
         call modin(filenam,pmodip)
         jdip=1
      endif
      if (ut1.ne.UT0.or.mth1.ne.mth0) then
         call sdec(mth1,ut1,sdelta,cdelta)
         UT0=ut1
         mth0=mth1
      endif
      if (flx1.ne.flx0) then
         R12=sqrt(167273.0D0+(flx1-63.7D0)*1123.6D0)-408.99D0
         flx0=flx1
      endif

      xlt=dmod(ut1+along/15.0D0+24.0D0,24.0D0)
      xMODIP=amodip(pmodip,alat,along)
      cchi=sin(alat*DR)*sdelta-cos(alat*DR)*cdelta*
     * cos(pi12*xlt)
      chi=atan2(sqrt(1.0D0-cchi*cchi),cchi)*RD
      call cciri(xMODIP,mth0,UT0,R12,alat,along,foF2,M3000)
      call ef1r(alat,mth1,flx1,chi,foF2, efoE,efoF1)
      call prepmdgr(mth0,R12,foF2,efoF1,efoE,M3000,hm,bb,aep)
      return
      end

      real*8 function  NeNeQ(h,hm,aep,bb)
      implicit real*8 (a-h,o-z)
      real*8 NeMdGR
      dimension hm(3),aep(3),bb(6)
 
      if (h.gt.hm(1)) then
         aNmax=NeMdGR(aep,hm,bb,hm(1))
         NeNeQ=topq(h,aNmax,hm(1),bb(6))
         return
      endif
      NeNeQ=NeMdGR(aep,hm,bb,h)
      return
      end

      real*8 function NeMdGR(aep,hm,bb,h)
      implicit real*8 (a-h,o-z)
      dimension aep(3),hm(3),bb(6),B(3)

      save

      parameter (f1=10.0D0, f2=1.0D0)
      parameter (h0=90.0D0, Hd=5.0D0)
      B(1)=bb(5)
      B(2)=bb(3)
      B(3)=bb(1)
      if (h.gt.hm(3)) B(3)=bb(2)
      if (h.gt.hm(2)) B(2)=bb(4)
         sum=0.0D0
         do jj=1,3
            arg0=(h-hm(jj))
            arg=arg0/B(jj)
            if (jj.gt.1) then
               d=abs(h-hm(1))
               arg=arg*exp(f1/(1.0D0+f2*d))
            endif
          if (h.lt.h0) arg=arg*(Hd+h0-h)/Hd
            if (abs(arg).gt.25.0D0) then
               s0=0.0D0
            else
               ee=exp(arg)
               s0=aep(jj)*ee/(1.0D0+ee)**2
            endif
            sum=sum+s0
         enddo
         NeMdGR=sum*1.0D11
         return
      end

      real*8 function topq(h,aNo,hmax,Ho)
      implicit real*8 (a-h,o-z)
      parameter (g=0.125D0,rfac=100.0D0)
         dh=h-hmax
         gh=g*dh
         z=dh/(Ho*(1.0D0+rfac*gh/(rfac*Ho+gh)))
         ee=fexp(z)
         if (ee.gt.1.0D11) then
            ep=4.0D0/ee
         else
            ep=4.0D0*ee/(1.0D0+ee)**2
         endif
         topq=aNo*ep
      return
      end

      subroutine prepmdgr(mth,R12,foF2,efoF1,efoE,M3000,hm,bb,aep)
      implicit real*8 (a-h,o-z)
      real*8 NmE,NmF1,NmF2,M3000
      dimension aep(3),hm(3),bb(6)
      data hmE /120.0D0/

      FNe(X)=0.124D0*X*X
      FEpst(X,Y,Z,W)=X*fexp((W-Y)/Z)/(1.0D0+fexp((W-Y)/Z))**2

      NmF2=FNe(foF2)
      NmF1=FNe(efoF1)
      NmE= FNe(efoE)
      hmF2=peakh(efoE,foF2,M3000)
      hmF1=0.5D0*(hmF2+hmE)
      hm(1)=hmF2
      hm(2)=hmF1
      hm(3)=hmE

      dNdHmx=-3.467D0+1.714D0*log(foF2)+2.02D0*log(M3000)
      dNdHmx=exp(dNdHmx)
      B2bot=0.385D2*NmF2/dNdHmx
      hF2F1E=hmF2-hmF1
      B1top=0.3D0*hF2F1E
      B1bot=0.5D0*hF2F1E
      Betop=B1bot
      if (Betop.lt.7.0D0) Betop=7.0D0
      Bebot=5.0D0
      aep(1)=4.0D0*NmF2
      aep(2)=4.0D0*NmF1
      aep(3)=4.0D0*NmE
      if (efoF1.le.0.5D0) then
         aep(2)=0.0D0
         aep(3)=4.0D0*(NmE-FEpst(aep(1),hmF2,B2bot,hmE))
      else
         do izml=1,5
            aep(2)=4.0D0*(NmF1-FEpst(aep(1),hmF2,B2bot,hmF1)
     &                      -FEpst(aep(3),hmE, Betop,hmF1))
            aep(2)=djoin(aep(2),0.8D0*NmF1,1.0D0,aep(2)
     &                                     -0.8D0*NmF1)
            aep(3)=4.0D0*(NmE -FEpst(aep(2),hmF1,B1bot,hmE)
     &                      -FEpst(aep(1),hmF2,B2bot,hmE))
         enddo
         if (aep(2).lt.0.0D0) then
            aep(2)=0.0D0
            aep(3)=4.0D0*(NmE-FEpst(aep(1),hmF2,B2bot,hmE))
         endif
      endif
      aep(3)=djoin(aep(3),0.5D-2,60.0D0,aep(3)-0.5D-2)
      b2k = 3.22D0-0.538D-1*foF2-0.664D-2*hmF2+0.113D0*hmF2/B2bot
     +    +0.257D-2*R12
      b2k=djoin(b2k,1.0D0,2.0D0,b2k-1.0D0)
      B2top=b2k*B2bot
      bb(1)=Bebot
      bb(2)=Betop
      bb(3)=B1bot
      bb(4)=B1top
      bb(5)=B2bot
      bb(6)=B2top
      return
      end

      subroutine ef1r(alat,mth,flx,chi,foF2,efoE,efoF1)
      implicit real*8 (a-h,o-z)
      parameter (DR=1.74532925199433D-2)
      parameter (chi0=86.23292796211615D0)

      chin=djoin(90.0D0-0.24D0*fexp(20.0D0-0.20D0*chi),chi,12.0D0,
     & chi-chi0)
      ee=fexp(0.3D0*alat)
      seas=1.0D0-2.0D0/(ee+1.0D0)
      if(mth.eq.1.or.mth.eq.2.or.mth.eq.11.or.mth.eq.12) 
     +  sfac=(1.112D0+0.019D0*seas)*flx**0.25D0
      if(mth.eq.3.or.mth.eq.4.or.mth.eq.9.or.mth.eq.10)
     +  sfac=1.112D0*flx**0.25D0
      if(mth.eq.5.or.mth.eq.6.or.mth.eq.7.or.mth.eq.8)
     +  sfac=(1.112D0-0.019D0*seas)*flx**0.25D0
      fa=sfac*fexp(log(cos(chin*DR))*0.3D0)
      efoE=sqrt(fa*fa+0.49D0)
      efoF1=djoin(1.4D0*efoE,0.0D0,1000.0D0,efoE-2.0D0)
      efoF1=djoin(0.0D0,efoF1,1000.0D0,efoE-efoF1)
      efoF1=djoin(efoF1,0.85D0*efoF1,60.0D0,0.85D0*foF2-efoF1)
      if (efoF1.lt.1.0D-6)efoF1=0.0D0
      return
      end

      subroutine cciri(xMODIP,mth,UT,R12,alat,along,foF2,M3000)
      implicit real*8 (a-h,o-z)
      real*8 M3000
      character*10 filena
      character*2 cm
      dimension FF0(988),xm0(441),F2(13,76,2),FM3(9,49,2)
      integer QM(7),QF(9)
      save
      data QF/11,11,8,4,1,0,0,0,0/,QM/6,7,5,2,1,0,0/
      data montha,monthb,Rga/13,14,-10.0D0/
      if (mth.ne.montha) then
         write(cm,'(I2.2)') mth+10
         filena='ccir'//cm//'.asc'
         open(77,file=filena,status='OLD')
         read(77,'(4e16.8)') F2,FM3
         close(77)
         montha=mth
      endif
      if (R12.ne.Rga.or.mth.ne.monthb) then
         RR2=R12/100.0D0
         RR1=1.0D0-RR2
         do i=1,76
         do j=1,13
            k=j+13*(i-1)
            FF0(k)=F2(j,i,1)*RR1+F2(j,i,2)*RR2
         enddo
         enddo
         do i=1,49
         do j=1,9
            k=j+9*(i-1)
            xm0(k)=FM3(j,i,1)*RR1+FM3(j,i,2)*RR2
         enddo
         enddo
         Rga=R12
         monthb=mth
      endif
      foF2= gamma1(xMODIP,alat,along,UT,6,QF,9,76,13,988,FF0)
      M3000=gamma1(xMODIP,alat,along,UT,4,QM,7,49, 9,441,xm0)
      return
      end

      real*8 function gamma1(xMODIP,alat,along,hour,iharm,nq,
     ,                          k1,m,mm,m3,sfe)
      implicit real*8 (a-h,o-z)
      real*8 c(12),s(12),coef(100),sum
      dimension nq(k1),xsinx(13),sfe(m3)
      logical numok
      parameter (DR=1.74532925199433D-2)

      hou=(15.0D0*hour-180.0D0)*DR
      s(1)=sin(hou)
      c(1)=cos(hou)
      do i=2,iharm
         c(i)=c(1)*c(i-1)-s(1)*s(i-1)
         s(i)=c(1)*s(i-1)+s(1)*c(i-1)
      enddo
      do i=1,m
         mi=(i-1)*mm
         coef(i)=sfe(mi+1)
         do j=1,iharm
            coef(i)=coef(i)+sfe(mi+2*j)*s(j)+sfe(mi+2*j+1)*c(j)
         enddo
      enddo
      sum=coef(1)
      ss=sin(xMODIP*DR)
      s3=ss
      xsinx(1)=1.0D0
      index=nq(1)
      do j=1,index
         numok=abs(ss).ge.1.0D-30
         if (numok) then
            sum=sum+coef(1+j)*ss
            xsinx(j+1)=ss
            ss=ss*s3
         else
            xsinx(j+1)=0.0D0
         endif
      enddo
      if (numok) then
         xsinx(nq(1)+2)=ss
      else
         xsinx(nq(1)+2)=0.0D0
      endif
      np=nq(1)+1
      ss=cos(alat*DR)
      s3=ss
      do j=2,k1
         s0=along*(j-1)*DR
         s1=cos(s0)
         s2=sin(s0)
         index=nq(j)+1
         do L=1,index
            np=np+1
            sum=sum+coef(np)*xsinx(L)*ss*s1
            np=np+1
            sum=sum+coef(np)*xsinx(L)*ss*s2
         enddo
         ss=ss*s3
      enddo
      gamma1=sum
      return
      end

      real*8 function peakh(foE,foF2,M3000)
      implicit real*8 (a-h,o-z)
      real*8 MF,M3000
      sqM=M3000*M3000
      MF=M3000*sqrt((0.0196D0*sqM+1.D0)/(1.2967D0*sqM-1.0D0))
      If(foE.ge.1.0D-30) then
         ratio=foF2/foE
         ratio=djoin(ratio,1.75D0,20.0D0,ratio-1.75D0)
         dM=0.253D0/(ratio-1.215D0)-0.012D0
      else
         dM=-0.012D0
      endif
      peakh=1490.0D0*MF/(M3000+dM)-176.0D0
      return
      end

      subroutine sdec(mth,UT,sdelta,cdelta)
      implicit real*8 (a-h,o-z)
      parameter (DR=1.74532925199433D-2)

      doy=mth*30.5D0-15.0D0
      t =doy + (18.0D0-UT)/24.0D0
      amrad=(0.9856D0*t - 3.289D0)*DR
      aLrad = amrad + (1.916D0*sin(amrad)+0.020D0*sin(2.0D0*amrad)+
     + 282.634D0)*DR
      sdelta=0.39782D0*sin(aLrad)
      cdelta=sqrt(1.0D0-sdelta*sdelta)
      return
      end

      subroutine modin(filenam,pmodip)
      implicit real*8 (a-h,o-z)
      character*15 filenam
      dimension pmodip(0:183,-1:182)
      character*100 txt
      parameter(latp=180,lngp=180,lathp=90,lnghp=90)

      open(31,status='OLD',file=filenam,form='FORMATTED')
      read(31,'(A)')txt
      do i=-lathp,lathp
         read(31,*) (pmodip(i+lathp+1,j+lnghp),j=-lnghp,lnghp)
      enddo
      close(31)
      do i=0,lngp
         pmodip(0,i)=pmodip(2,mod((i+lnghp),lngp))
         pmodip(latp+2,i)=pmodip(latp,mod((i+lnghp),lngp))
         pmodip(latp+3,i)=pmodip(latp-1,mod((i+lnghp),lngp))
      enddo
      do i=0,latp+3
         pmodip(i,-1)=pmodip(i,lngp-1)
         pmodip(i,lngp+1)=pmodip(i,1)
      enddo
      return
      end


      real*8 function amodip(pmodip,alat,along)
      implicit real*8 (a-h,o-z)
      dimension pmodip(0:183,-1:182)
      dimension z(4),z1(4)
      parameter(lngp=180,dlatp=1.0D0,dlngp=2.0D0)

      along1=dmod(along+360.0D0,360.0D0)
      dlng1=(along1+180.0D0)/dlngp
      dlng1=dlng1-dint(dlng1)
      j1=idint((along1+180.0D0)/dlngp)-2
      if (j1.lt.0) j1=j1+lngp
      if (j1.gt.lngp-3) j1=j1-lngp
      a=(alat+90.0D0)/dlatp+1.0D0
      i=idint(a)-2
      a=a-dfloat(i+2)
      do k = 1,4
         do j=1,4
         z1(j)=pmodip(i+j,j1+k)
         enddo
         z(k)=finter3(z1,a)
      enddo
      amodip=finter3(z,dlng1)
      return
      end

      real*8 function fexp(a)
      real*8 a
      if(a.gt.80.0D0) then
         fexp=5.5406D34
         return
      endif
      if(a.lt.-80.0D0) then
         fexp=1.8049D-35
         return
      endif
      fexp=exp(a)
      return
      end

      real*8 function djoin(f1,f2,alpha,x)
      real*8 f1,f2,alpha,x,ee,fexp
      ee=fexp(alpha*x)
      djoin=(f1*ee+f2)/(ee+1.0D0)
      return
      end

      real*8 function finter3(z,x)
      implicit real*8 (a-h,o-z)
      dimension z(4),a(0:3)

      dx=x+x-1.0D0
      if (abs(dx+1.0D0).lt.1.0D-10) then
         finter3=z(2)
      else
         g1=(z(3)+z(2))
         g2=(z(3)-z(2))
         g3=(z(4)+z(1))
         g4=(z(4)-z(1))/3.0D0
         a(0)=(9.0D0*g1-g3)
         a(1)=(9.0D0*g2-g4)
         a(2)=(g3-g1)
         a(3)=(g4-g2)
         zi=0.0D0
         do j=3,0,-1
            zi=zi*dx+a(j)
         enddo
      finter3=zi/16.0D0
      endif
      return
      end
