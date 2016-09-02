C=========================================================================
C     NeQVal for NeQuick2 P.531 
C     Release date 22 Oct 2011
C     (Internal developers reference 2.1.0)
C ========================================================================
C DISCLAIMER:
C	This software is meant for scientific use only.
C         Please acknowledge the Aeronomy and Radiopropagation Laboratory
C     	of the Abdus Salam International Centre for Theoretical Physics
C     	Trieste, Italy
C
C     This software is provided to The User “as is”. ITU-R and the authors assume 
C	no responsibility whatsoever for its use by other parties, and makes no 
C	guarantees, expressed or implied, about the quality, reliability, or any 
C 	other characteristic of this software. Under no circumstances shall ITU, 
C	the author(s) or contributor(s) be liable for damages resulting directly or 
C	indirectly from the use, misuse or inability to use this software.
C	
C    The user may not modify this software in content and then present it or 
C	results derived from it as ITU-R or author(s) material. 
C =========================================================================
C
C usage: NeQVal F10.7 InputFile (> OutputFile) 
C calculates Slant Total Electron Content for a given Solar Flux at 
C 10.7 cm (F10.7) all lines in InputFile, where 
C each line defines the epoch and the coordinates of begin and end 
C points, the year. The results is provided in the Standard Output 
C with a line including the same line of the input file plus an addtional 
C column with the Slant Total Electron Content along the ray.
C
C    Output file: selected by user
C
C =========================================================================
C
C Output (function value) is electron density in TECU 
C   (Total Electron Content Units)
C 
C =========================================================================
C     compilation: link object code with object code of NeQuick_2
C	information on compilation available in COMPILING.txt
C= ========================================================================
C
C   subroutines
C     integration routines based on slQu.for source from Dr. R. Leitinger.
C
C========================================================================
C     Driver progamme based on slQu_2 driver (Dr. Reinhart Leitinger et al)
C 
C     Adaptations from slQu_2 to NeQVal: 
C	Raul Orus-Perez, Roberto Prieto-Cerdeira
C	European Space Agency
C	Noordwijk, The Netherlands
C ========================================================================
 
      program NeQVal


      implicit double precision (a-h,o-z)
      parameter (Re=6371.2D0)
      character filein*256,cflux*5
      character fihelp*256


      ninput=iargc()
ccc Checking for the correct number of input arguments

      if (ninput.ne.2) then
ccc Displaying help
        call getarg(0,fihelp)

        write(*,*)"Usage: ",
     2fihelp(1:index(fihelp,'NeQuick_test_v0')+15)
     3," <File Input> flux"
        write(*,*)"Ex: ",fihelp(1:index(fihelp,'NeQuick_test_v0')+15)
     4,"/home/ITU/test/vill_ne_in.dat 70"
        stop
      endif
ccc Reading the input arguments
      call getarg(1,filein)
      call getarg(2,cflux)
      read(cflux,*)flx

      if (flx.gt.193.or.flx.lt.63) goto 870      

      open(10,file=filein,status="old",err=890)

10    read(10,*,end=99,err=8901)iyr,mth,UT,alng1,ph1,h1,alng2,ph2,h2

      h1=h1/1000.0d0
      h2=h2/1000.0d0


      call NeQuistec(ph1,alng1,h1,ph2,alng2,h2,mth,UT,flx,stec)
    
      write(*,123)iyr,mth,UT,alng1,ph1,h1,alng2,ph2,h2,stec
123   format (i4,1X,i2,1X,f10.5,1X,8(f16.5,1X))

      goto 10
99    continue 

      stop
ccc Managing some errors
870   stop 'Flux out of boundary; 63 < flx < 193 '
890   stop 'Input file is missing'
8901  stop 'Input data is corrupted'
      end

      subroutine NeQuistec(ph1,alng1,h1,ph2,alng2,h2
     1,mth1,UT1,flx1,stec)


      implicit real*8 (a-h,o-z)
C      real*8 NeQuick
C      real*8 NeNeQ
      dimension hm(3),aep(3),bb(6)
      common /timin/flx,UT,mth
      common /nepar/hm,aep,bb
      parameter (Re=6371.2D0)
      parameter (RD=5.729577951308232D1)

      flx=flx1
      mth=mth1
      UT=UT1
 
      call rays(r1,h1,ph1,alng1,r2,h2,ph2,alng2,zeta,
     & pp,Re,sa,ca,sb,cb,ssig,csig,along1)
      s1=sqrt(r1*r1-pp*pp)
      s2=sqrt(r2*r2-pp*pp)

      R12=sqrt(167273.0D0+(flx-63.7D0)*1123.6D0)-408.99D0 

      h0=0.0D0
      if (h1.gt.h0) h0=h1
      r0=Re+h0
      s0=sqrt(r0*r0-pp*pp)
C  use of NeQuick necessary to condition the model to   mth,flx,UT
      call prepNeQ(ph1,alng1,mth,UT,flx,hm,aep,bb)

      if (pp.lt.0.1D0) then
C  integration along vertical profile
        alat=ph1
        along=alng1
        if (h2.le.1000.0D0) then
          tec1=gintv(h0,h2,1.0D-3)
        else
          h1a=1000.0D0
          if (h2.le.2000.0D0) then
            if (h1.ge.1000.0D0) then
              tec1=gintv(h1,h2,1.0D-3)
            else
              tec1=gintv(h0, h1a,1.0D-3)
              tec2=gintv(h1a,h2, 1.0D-2)
              tec4=tec1+tec2
            endif
          else
            if (h1.ge.2000.0D0) then
              tec1=gintv(h1,h2,1.0D-3)
            else
              h1b=2000.0D0
              if (h1.ge.1000.0D0) then
                tec1=gintv(h1,h1b,1.0D-3)
                tec2=gintv(h1b,h2,1.0D-3)
                tec4=tec1+tec2
              else
                tec1=gintv(h0, h1a,1.0D-3)
                tec2=gintv(h1a,h1b,1.0D-2)
                tec3=gintv(h1b,h2, 1.0D-2)
                tec4=tec1+tec2+tec3
              endif
            endif
          endif
        endif
      else
C  integration along slant profile
        if (h2.le.1000.0D0) then
           tec1=gint(s0,s2,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
        else
          h1a=1000.0D0
          r1a=h1a+Re
          s1a=sqrt(r1a*r1a-pp*pp)
          if (h2.le.2000.0D0) then
            if (h1.ge.1000.0D0) then
              tec1=gint(s1,s2,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
            else
              s2=sqrt(r2*r2-pp*pp)
              tec1=gint(s0, s1a,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
              tec2=gint(s1a,s2, 1.0D-2,  pp,Re,sa,ca,ssig,csig,along1)
              tec4=tec1+tec2
            endif
          else
            if (h1.ge.2000.0D0) then
              tec1=gint(s1,s2,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
            else
              h1b=2000.0D0
              r1b=h1b+Re
              s1b=sqrt(r1b*r1b-pp*pp)
              if (h1.ge.1000.0D0) then
                tec1=gint(s1,s1b,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
                tec2=gint(s1b,s2,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
                tec4=tec1+tec2
              else
                tec1=gint(s0, s1a,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
                tec2=gint(s1a,s1b,1.0D-2,  pp,Re,sa,ca,ssig,csig,along1)
                tec3=gint(s1b,s2, 1.0D-2,  pp,Re,sa,ca,ssig,csig,along1)
                tec4=tec1+tec2+tec3
              endif
            endif
          endif
        endif
      endif

      stec = tec4/1.0D13
      return
      end

      real*8 function gint (g1,g2,eps,pp,Re,s1,c1,ssig,csig,along1)
C  special; integrates ELD
      implicit real*8 (a-h,o-z)
         n = 8
    1    h  =  (g2-g1) / dfloat(n)
         hh = 0.5D0*h
         g  = h*0.5773502691896D0
         y  = g1 + (h-g)*0.5D0
         gint2  =  eld(y,pp,Re,s1,c1,ssig,csig,along1)+
     +             eld(y+g,pp,Re,s1,c1,ssig,csig,along1)
         do m = 1,n-1
            gint2 = gint2 + eld(y+h+g,pp,Re,s1,c1,ssig,csig,along1)+
     +                      eld(y+h,pp,Re,s1,c1,ssig,csig,along1)
            y = y + h
         enddo
         gint2 = gint2*hh
         if (n.eq.8.or.abs(gint1-gint2).gt.eps*abs(gint1)) then
            n = n*2
            gint1 = gint2
            if (n.lt.1024) goto 1
         endif
         gint = gint2+(gint2-gint1)/15.0D0
      return
      end

      real*8 function gintv (g1,g2,eps)
C  special; integrates NeNeQ for vertical profile
      implicit real*8 (a-h,o-z)
      real*8 NeNeQ
      dimension hm(3),aep(3),bb(6)
      common /nepar/hm,aep,bb
         n = 8
    1    h  =  (g2-g1) / dfloat(n)
         hh = 0.5D0*h
         g  = h*0.5773502691896D0
         y  = g1 + (h-g)*0.5D0
         gint2  =  NeNeQ(y,hm,aep,bb)+
     +             NeNeQ(y+g,hm,aep,bb)
         do m = 1,n-1
            gint2 = gint2 + NeNeQ(y+h+g,hm,aep,bb)+
     +                      NeNeQ(y+h,hm,aep,bb)
            y = y + h
         enddo
         gint2 = gint2*hh
         if (n.eq.8.or.abs(gint1-gint2).gt.eps*abs(gint1)) then
            n = n*2
            gint1 = gint2
            if (n.lt.1024) goto 1
         endif
         gintv = gint2+(gint2-gint1)/15.0D0
      return
      end

      subroutine rays(r1,h1,ph1,alng1,r2,h2,ph2,alng2,zeta,
     & pp,Re,sp,cp,s2,c2,ssig,csig,along1)
      implicit real*8 (a-h,o-z)
      akappa=Re/(Re+400.0D0)
   10 continue
c       write(6,*)'INPUT: ',
c     &'Ray endpoint 1: latitude (deg N), longitude (deg E), height (km)'
c      read(5,*)ph1,alng1,h1
c      write(6,*)'INPUT: ',
c     &'Ray endpoint 2: latitude (deg N), longitude (deg E), height (km)'
c      read(5,*)ph2,alng2,h2
C
C put "vertical" for "near vertical"
      if (abs(ph2-ph1).lt.1.0D-5.and.abs(alng2-alng1).lt.1.0D-5) then
         ph2=ph2
         alng2=alng1
      endif
      r1=Re+h1
      r2=Re+h2
C
C provides coordinates (pp,php,alamp) of ray perigee,
C zenith angle (zeta) for the ray endpoint 1,
C slant to vertical projection factor (cchi) for the
C conditions given by akappa (<1)
      call naut(r1,r2,ph1,ph2,alng1,alng2,akappa,
     &  pp,php,alamp,zeta,cchi)
      if (abs(zeta).gt.90.0.and.pp.lt.Re) then
         write(*,*) ' ray cuts surface of Earth'
         write(*,*) '      or endpoint 2 lower than endpoint 1.'
         write(*,*) ' Repeat input'
         goto 10
      endif
C parameters of ray properties for use in other modules (formerly in
C    a common block)
C    (origin for the s coordinate is the ray perigee with
C       lat-long coordinates php and alamp and the radius pp)
C    sp,cp:     sine and cosine of latitude point p (= ray perigee)
C    s2,c2:     sine and cosine of latitude point 2
C    ssig,csig: sine and cosine of ray azimuth
C               (point 2 seen from point 1)
      if (pp.ge.0.1D0) then
       call gcirc(php,ph2,alamp,alng2,sp,cp,s2,c2,ssig,csig,psi)
      endif
C    along1: longitude of point 1 (= ray perigee)
      along1=alamp
      return
      end

      subroutine gcirc(alat1,alat2,along1,along2,s1,c1,s2,c2,ssig,
     & csig,psi)
C
C  calculates great circle path properties
C  input parameters are the endpoints coordinates
C    (alat1,along1), (alat2,along2) [deg.]
C  output parameters are sine and cosine of the endpoint latitudes
C    s1,c1 and s2,c2, resp.,
C    sine and cosine of the azimuth of endpoint 2 seen from endpoint 1
C    ssig, csig (N over E to S)
C    and the great circle distance psi [deg.]
C
      implicit real*8 (a-h,o-z)
      parameter (DR=1.74532925199433D-2)
      rlat1=alat1*DR
      rlat2=alat2*DR
      dlong=(along2-along1)*DR
      s1=sin(rlat1)
      s2=sin(rlat2)
      c1=cos(rlat1)
      c2=cos(rlat2)
      sd=sin(dlong)
      cd=cos(dlong)
      if (abs(abs(alat1)-90.0D0).lt.1.0D-10) then
         psi=abs(alat2-alat1)
         ssig=0.0D0
         if (alat1.gt.0.0D0) then
            csig=-1.0D0
         else
            csig=1.0D0
         endif
      else
         cpsi=s1*s2+c1*c2*cd
         spsi=sqrt(1.0D0-cpsi*cpsi)
         ssig=c2*sd/spsi
         csig=(s2-s1*cpsi)/c1/spsi
         psi=atan2(spsi,cpsi)/DR
      endif
      return
      end

      subroutine naut(r1,r2,ph1,ph2,alng1,alng2,akappa,
     & pp,php,alamp,zeta,cchi)
C
C *** new version:
C  additional input parameter akappa, add. output parameter cchi ***
C
C  calculates position of ray perigee, zenith angle of ray at
C  lower endpoint (endpoint 1), and slant to vertical projection
C  factor cos(chi)
C
C  input parameters are the endpoint coordinates
C    (r1,ph1,alng1) (r2,ph2,alng2)
C    (distance from Earth centre [km], latitude, longitude [deg.])
C    and the ratio akappa=r1/(r1+hi), hi being the "mean ionospheric
C    height" if endpoint 1 is at the surface of the Earth
C  output parameters are the coordinates of the ray perigee
C    (point on the ray [straight line through the endpoints] closest
C     to the centre of the Earth)
C    (pp,php,alamp),
C    the zenith angle of endpoint 2 seen from entpoint 1, zeta [deg.]
C    and the projection factor cchi=cos(chi), chi being the zenith
C    angle in the point with height hi above endpoint 1.
C
      implicit real*8 (a-h,o-z)
      parameter (DR=1.74532925199433D-2,RD=5.729577951308232D1)
      parameter ( pi=3.141592653589793D0)
      if (abs(ph1-ph2).lt.1.0D-5.and.abs(alng1-alng2).lt.1.0D-5) then
C vertical profile
         pp=0.0D0
         php=ph1
         alamp=alng1
         zeta=0.0D0
         cchi=1.0D0
      else
C slant profile
         sph1=sin(ph1*DR)
         cph1=cos(ph1*DR)
         sph2=sin(ph2*DR)
         cph2=cos(ph2*DR)
         cdl12=cos((alng2-alng1)*DR)
         sdl12=sin((alng2-alng1)*DR)
         cdel=sph1*sph2+cph1*cph2*cdl12
         sdel=sqrt(1.0D0-cdel*cdel)
         zeta=atan2(sdel,cdel-r1/r2)
         ssigp=sdl12*cph2/sdel
         csigp=(sph2-cdel*sph1)/sdel/cph1
         delp=-zeta+pi/2.0D0
         sdelp=sin(delp)
         cdelp=cos(delp)
         sphp=sph1*cdelp-cph1*sdelp*csigp
         cphp=sqrt(1.0D0-sphp*sphp)
         php=atan2(sphp,cphp)*RD
         slamp=-ssigp*sdelp/cphp
         clamp=(cdelp-sph1*sphp)/cph1/cphp
         alamp=atan2(slamp,clamp)*RD+alng1
         szeta=sin(zeta)
         pp=r1*szeta
         zeta=zeta*RD
         schi=akappa*szeta
         cchi=sqrt(1.0D0-schi*schi)
      endif
      return
      end

      subroutine geogra(s, pp,Re,s1,c1,ssig,csig,along1, h,alat,along)
C
C  calculates height, latitude and longitude along the given "ray"
C  for given distance from the perigee of the ray
C  input parameter: distance from the perigee of the ray, s (km)
C  output parameters:
C    height h (km), latitude alat (deg. N), longitude along (deg. E)
C
C  the properties of the ray pp,s1,c1,ssig,csig,along1
C      and the Earth radius Re are input parameters
C
      implicit real*8 (a-h,o-z)
      parameter (RD=5.729577951308232D1)
      tdel=s/pp
      cdel=1.0D0/sqrt(1.0D0+tdel*tdel)
      sdel=tdel*cdel
      arg=s1*cdel+c1*sdel*csig
      alat=atan2(arg,sqrt(1.0D0-arg*arg))*RD
      clong=atan2(sdel*ssig*c1,cdel-s1*arg)*RD
      along=clong+along1
      h=sqrt(s*s+pp*pp)-Re
      return
      end

      real*8 function eld(s,pp,Re,s1,c1,ssig,csig,along1)
C
C  gives electron density as a function of the coordinate s
C    (distance from the perigee of the ray)
C
      implicit real*8 (a-h,o-z)
      real*8 NeQuick
      common /timin/flx,UT,mth
      parameter (RD=5.729577951308232D1)
      tdel=s/pp
      cdel=1.0D0/sqrt(1.0D0+tdel*tdel)
      sdel=tdel*cdel
      arg=s1*cdel+c1*sdel*csig
      alat=atan2(arg,sqrt(1.0D0-arg*arg))*RD
      clong=atan2(sdel*ssig*c1,cdel-s1*arg)*RD
      along=clong+along1
      h=sqrt(s*s+pp*pp)-Re
      eld=NeQuick(h,alat,along,mth,flx,UT)
      return
      end
