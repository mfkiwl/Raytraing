C     ========================================================================
C     slQu for NeQuick2 P.531
C     Release date 22 Oct 2011
C     (Internal developers reference 2.1.0)
C     ========================================================================
C     DISCLAIMER:
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
C     =========================================================================
C     
C     Main program to calculate electron content from NeQuick2 P.531
C     for arbitrarily chosen rays which do not cut the surface of
C     the Earth between the given endpoints.
C     
C     Output file: slQu.dat
C     
C     This file (main program plus auxiliary modules listed below)
C     is written in FORTRAN 77.
C     Implicit double precision for all real variables and
C     functions (real*8) is used.
C     It consists of 558 FORTRAN statement lines
C     (comment and empty lines excluded).
C     
C     ========================================================================
C     compilation: link object code with object code of NeQuick_2
C	information on compilation available in COMPILING.txt
C     ========================================================================
C     
C     ray conventions:
C     spherical Earth
C     straight line "rays"
C     coordinate s [km] along the ray, origin in ray perigee, point of
C     ray closest to the centre of the Earth.
C     It has radius pp [km].
C     Latitude and longitude of ray perigee: php, alamp
C     input of a ray: subroutine RAYS uses the endpoint coordinates
C     (height in km, latitude in deg. N, longitude in deg. E)
C     to define a ray (1 for lower, 2 for higher endpoint)
C     
C     ========================================================================
C     additional data file (optional): R12.dat
C     each line contains the year (format I5) and the 12 R12 data for this year
C     (format 12 F6.1), e.g.
C     2000  112.9 116.8 119.9 120.8 119.0 118.7 119.8 118.6 116.3 114.5 112.7 112.0
C     2001  108.7 104.0 104.8 107.6 108.6 109.8 111.7 113.6 114.1 114.0 115.5 114.9
C     2001  108.7 104.0 104.8 107.5 108.6 109.8 111.7 113.6 114.1 114.0 115.5 114.6
C     2002  113.5 114.6 113.3 110.5 108.8 106.2 102.7  98.7  94.6  90.5  85.2  82.0
C     2003   80.8  78.3  74.0  70.1  67.6  65.0  61.8  60.0  59.5  58.2  56.7  54.8
C     2004   52.0  49.3  47.1  45.5  43.8  41.6  40.2  39.2  37.5  35.9  35.3  35.2
C     2005   34.6  33.9  33.5  31.6  28.9  28.8  29.1  27.4  25.8  25.5  24.9  23.0
C     ========================================================================
C     
C     This source code contains the following modules:
C     
C     (1) numerical integration along a slant ray
C     
C     real*8 function gint (g1,g2,eps,pp,Re,s1,c1,ssig,csig,along1)
C     numerical integration (specialized, no external function)
C     
C     (2) numerical integration along a vertical ray
C     
C     real*8 function gintv (g1,g2,eps)
C     numerical integration (specialized, no external function)
C     
C     (3) separation in functional units
C     subroutine rays(r1,h1,ph1,alng1,r2,h2,ph2,alng2,zeta,
C     & pp,Re,sp,cp,s2,c2,ssig,csig,along1)
C     set and check ray endpoints, calculate geometric parameters for ray
C     
C     subroutine dat_t_sa(iyr,mth,nday,tt,R12,flx)
C     set date and solar activity
C     
C     (4) package with auxiliary subroutines and functions
C     subroutine gcirc(alat1,alat2,along1,along2,s1,c1,s2,c2,ssig,
C     & csig,psi)
C     calculates great circle path properties
C     
C     subroutine naut(r1,r2,ph1,ph2,alng1,alng2,akappa,
C     & pp,php,alamp,zeta,cchi)
C     calculates position of ray perigee, zenith angle of ray at
C     lower endpoint, and slant to vertical projection factor cos(chi)
C     
C     subroutine geogra(s,pp,Re,s1,c1,ssig,csig,along1,
C     &   h,alat,along)
C     calculates height, latitude and longitude along the given "ray"
C     for given distance from the perigee of the ray, s
C     
C     real*8 function eld(s,pp,Re,s1,c1,ssig,csig,along1)
C     gives electron density as a function of the coordinate s
C     
C     gint, rays, geogra, eld: input parameters pp,Re,s1,c1,ssig,csig,along1
C     they provide the ray properties
C     pp:     radius of perigee
C     Re:     Earth radius
C     s1:     sine   of latitude of lower end point
C     c1:     cosine of latitude of lower end point
C     ssig:   sine   of azimuth of ray at lower end point
C     csig:   cosine of azimuth of ray at lower end point
C     along1: longitude of lower end point
C     
C     ========================================================================
C     
C     Author of the original driver program:
C     Dr. Reinhart Leitinger
C     Responsible authors for this package:
C     Bruno Nava and Yenca Migoya Orue'
C     the Abdus Salam International Centre for Theoretical Physics
C     strada costiera 11
C     34014 Trieste (TS)
C     Italy
C     e-mail: bnava@ictp.it
C     yenca@ictp.it
C     
C     with contributions by Johanna Wagner (IGAM) and Pierdavide Coisson.
C     
C     ========================================================================
C     
C     Changes:
C     14.11.2007 all constants in double precision R12.dat includes 2006
C     11.02.2008 output file includes always the units for TEC
C     08.06.2010 R12.dat extended to include 2009 data
C     21.06.2010 ctrl of user input of flx and R12

      program slQu
      implicit real*8 (a-h,o-z)
C     real*8 NeQuick
      real*8 NeNeQ
      character*80 filen1
      character*1 yn
      dimension hm(3),aep(3),bb(6)
      common /timin/flx,UT,mth
      common /nepar/hm,aep,bb
      parameter (Re=6371.2D0)
      parameter (RD=5.729577951308232D1)

      write(6,*)
      write(6,*)'          *****************************************'
      write(6,*)'          *             NeQuick2 P.531            *'
      write(6,*)'          *   slant profile and electron content  *'
      write(6,*)'          *                                       *'
      write(6,*)'          * This software is meant for scientific *'
      write(6,*)'          * use only.                             *'
      write(6,*)'          * Conditions of use in Readme.txt       *'
      write(6,*)'          *                                       *'
      write(6,*)'          *****************************************'
      write(6,*)
      write(6,*)
     &     'Electron density is calculated along straight line rays'
      write(6,*)
     &     '   from a lower endpoint (1) to a higher one (2).'
      write(6,*)

      filen1='slQu.dat'
      open(16,file=filen1)

      call rays(r1,h1,ph1,alng1,r2,h2,ph2,alng2,zeta,
     &     pp,Re,sa,ca,sb,cb,ssig,csig,along1)
      s1=sqrt(r1*r1-pp*pp)
      s2=sqrt(r2*r2-pp*pp)
      write(16,'(A/2F8.2,F9.2)')
     &     'Ray endpoint 1: lat. (deg. N), long. (deg. E), height (km)',
     &     ph1,alng1,h1
      write(16,'(A/2F8.2,F9.2)')
     &     'Ray endpoint 2: lat. (deg. N), long. (deg. E), height (km)',
     &     ph2,alng2,h2
      if (pp.ge.0.1D0) then
         call gcirc(ph1,ph2,alng1,alng2,sp1,cp1,sp2,cp2,saz,caz,psi)
         write(16,'(2A/2F8.2)')
     &        'zenith angle (deg.) and azimuth (N over E to S, deg.)',
     &        ' of ray at endpoint 1 ',zeta,atan2(saz,caz)*RD
      endif

      call dat_t_sa(4,iyr,mth,nday,ut,R12,flx)

      write(16,'(A,I4,1H,2(F6.1,1H,),I3,1H,,F5.1)')
     +     'Year, S10.7, R12, month, UT: ',iyr,flx,R12,mth,ut
      write(16,'(/A/2A)')
     &     'Electron contents along ray.',
     &     '  (h1-h2) means from point in ',
     &     'height h1 to point in height h2 (heights in km)'

      write(6,*)'List electron density profile along ray (y/n)?'
      read(5,'(A)')yn
      h0=0.0D0
      if (h1.gt.h0) h0=h1
      r0=Re+h0
      s0=sqrt(r0*r0-pp*pp)
C     use of NeQuick necessary to condition the model to  mth,flx,UT
      call prepNeQ(ph1,alng1,mth,UT,flx,hm,aep,bb)
      if (yn.eq.'Y'.or.yn.eq.'y') then
         if (pp.ge.0.1D0) write(16,*) 's: coordinate along ray'
         write(16,*)
     &        'r: radius (distance from center of Earth)'
         write(16,*)
         if (pp.lt.0.1D0) then
            write(16,*)'     r    height  lat   long   el.density'
            write(16,*)'    km      km   deg N  deg E    m^-3'
         else
            write(16,*)'   s       r    height  lat   long   el.density'
            write(16,*)'  km      km      km   deg N  deg E    m^-3'
         endif
         dh=10.0D0
         if (h1.ge.500.0D0) dh= 50.0D0
         if (h1.ge.2000.0D0)dh=250.0D0
         h=h1-dh
 10      h=h+dh
         r=h+Re
         if (pp.lt.0.1D0) then
            aNe=NeNeQ(h,hm,aep,bb)
            if (aNe.lt.1000.0D0) aNe=0.0D0
            write(16,'(2F8.1,2F8.2,E13.6)') r,h,ph1,alng1,aNe
         else
            s=sqrt(r*r-pp*pp)
C     (call of geogra needed only to provide  alat,along  for output)
            call geogra (s,  pp,Re,sa,ca,ssig,csig,along1,h,alat,along)
            aNe=eld(s, pp,Re,sa,ca,ssig,csig,along1)
            if (aNe.lt.1000.0D0) aNe=0.0D0
            alongw=dmod(dmod(along+360.0D0,360.0D0)+360.0D0,360.0D0)
            write(16,'(3F8.1,2F8.2,E13.6)') s-s1,r,h,alat,alongw,aNe
         endif
         if (nint(h).eq.500)  dh=50.0D0
         if (nint(h).eq.2000) dh=250.0D0
         if (h+dh.le.h2) goto 10
         if (h+0.01D0.lt.h2) then
            r=h2+Re
            s=sqrt(r*r-pp*pp)
            if (pp.lt.0.1D0) then
               aNe=NeNeQ(h2,hm,aep,bb)
               if (aNe.lt.1000.0D0) aNe=0.0D0
               alongw=dmod(dmod(along+360.0D0,360.0D0)+360.0D0,360.0D0)
               write(16,'(2F8.1,2F8.2,E13.6)') r,h,alat,alongw,aNe
            else
C     (call of geogra needed only to provide  alat,along for output)
               call geogra (s, pp,Re,sa,ca,ssig,csig,along1,
     &              h,alat,along)
               aNe=eld(s, pp,Re,sa,ca,ssig,csig,along1)
               if (aNe.lt.1000.0D0) aNe=0.0D0
               alongw=dmod(dmod(along+360.0D0,360.0D0)+360.0D0,360.0D0)
               write(16,'(3F8.1,2F8.2,E13.6)') s-s1,r,h2,alat,alongw,aNe
            endif
         endif
      endif

      if (pp.lt.0.1D0) then
C     integration along vertical profile
         alat=ph1
         along=alng1
         if (h2.le.1000.0D0) then
            tec1=gintv(h0,h2,1.0D-3)
            write(16,'(A,I4,1H-,I4,1H))')
     &           'Electron content (',nint(h0),nint(h2)
            write(16,'(16X,F12.2,A)')tec1/1.0D12,'   x10^15 m^-2'
         else
            h1a=1000.0D0
            if (h2.le.2000.0D0) then
               if (h1.ge.1000.0D0) then
                  tec1=gintv(h1,h2,1.0D-3)
                  write(16,'(A,I4,1H-,I4,1H) )')
     &                 'Electron content (',nint(h1),nint(h2)
                  write(16,'(16X,F12.2,A)')tec1/1.0D12,'   x10^15 m^-2'
               else
                  tec1=gintv(h0, h1a,1.0D-3)
                  tec2=gintv(h1a,h2, 1.0D-2)
                  tec4=tec1+tec2
                  write(16,'(A,2(I4,1H-,I4,3H),(),I4,1H-,I4,1H))')
     &                 'Electron contents (',nint(h0),nint(h1a),
     &                 nint(h1a),nint(h2),
     &                 nint(h0), nint(h2)
                  write(16,'(16X,3F12.2,A)')
     &                 tec1/1.0D12,tec2/1.0D12,tec4/1.0D12,   
     &                 'x10^15 m^-2'
               endif
            else
               if (h1.ge.2000.0D0) then
                  tec1=gintv(h1,h2,1.0D-3)
                  write(16,'(A,I4,1H-,I5,1H) )')
     &                 'Electron content (',nint(h1),nint(h2)
                  write(16,'(16X,F12.2,A)')tec1/1.0D12,'   x10^15 m^-2'
               else
                  h1b=2000.0D0
                  if (h1.ge.1000.0D0) then
                     tec1=gintv(h1,h1b,1.0D-3)
                     tec2=gintv(h1b,h2,1.0D-3)
                     tec4=tec1+tec2
                     write(16,'(A,I4,1H-,I4,3H),(,I4,1H-,I5,3H),(,
     &I4,1H-,I5,1H))')
     &                    'Electron contents (',nint(h1),nint(h1b),
     &                    nint(h1b),nint(h2),
     &                    nint(h1), nint(h2)
                     write(16,'(16X,3F12.2,A)')
     &                    tec1/1.0D12,tec2/1.0D12,tec4/1.0D12,
     &                    '   x10^15 m^-2'
                  else
                     tec1=gintv(h0, h1a,1.0D-3)
                     tec2=gintv(h1a,h1b,1.0D-2)
                     tec3=gintv(h1b,h2, 1.0D-2)
                     tec4=tec1+tec2+tec3
                     write(16,'(A,2(I4,1H-,I4,3H),(),I4,1H-,I5,3H),(, 
     &I4,1H-,I5,1H))')
     &                    'Electron contents (',nint(h0),nint(h1a),
     &                    nint(h1a),nint(h1b),
     &                    nint(h1b),nint(h2),
     &                    nint(h0), nint(h2)
                     write(16,'(16X,4F12.2,A)')
     &                    tec1/1.0D12,tec2/1.0D12,tec3/1.0D12,
     &                    tec4/1.0D12,'   x10^15 m^-2'
                  endif
               endif
            endif
         endif
      else
C     integration along slant profile
         if (h2.le.1000.0D0) then
            tec1=gint(s0,s2,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
            write(16,'(A,I4,1H-,I4,1H))')
     &           'Electron content (',nint(h0),nint(h2)
            write(16,'(16X,F12.2,A)')tec1/1.0D12,'   x10^15 m^-2'
         else
            h1a=1000.0D0
            r1a=h1a+Re
            s1a=sqrt(r1a*r1a-pp*pp)
            if (h2.le.2000.0D0) then
               if (h1.ge.1000.0D0) then
                  tec1=gint(s1,s2,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
                  write(16,'(A,I4,1H-,I4,1H) )')
     &                 'Electron content (',nint(h1),nint(h2)
                  write(16,'(16X,F12.2,A)')tec1/1.0D12,'   x10^15 m^-2'
               else
                  s2=sqrt(r2*r2-pp*pp)
                  tec1=gint(s0, s1a,1.0D-3,pp,Re,sa,ca,ssig,csig,along1)
                  tec2=gint(s1a,s2, 1.0D-2,pp,Re,sa,ca,ssig,csig,along1)
                  tec4=tec1+tec2
                  write(16,'(A,2(I4,1H-,I4,3H),(),I4,1H-,I4,1H))')
     &                 'Electron contents (',nint(h0),nint(h1a),
     &                 nint(h1a),nint(h2),
     &                 nint(h0), nint(h2)
                  write(16,'(16X,3F12.2,A)')
     &                 tec1/1.0D12,tec2/1.0D12,tec4/1.0D12,
     &                 '   x10^15 m^-2'
               endif
            else
               if (h1.ge.2000.0D0) then
                  tec1=gint(s1,s2,1.0D-3,  pp,Re,sa,ca,ssig,csig,along1)
                  write(16,'(A,I4,1H-,I5,1H) )')
     &                 'Electron content (',nint(h1),nint(h2)
                  write(16,'(16X,F12.2,A)')tec1/1.0D12,'   x10^15 m^-2'
               else
                  h1b=2000.0D0
                  r1b=h1b+Re
                  s1b=sqrt(r1b*r1b-pp*pp)
                  if (h1.ge.1000.0D0) then
                     tec1=gint(s1,s1b,1.0D-3,  
     &                    pp,Re,sa,ca,ssig,csig,along1)
                     tec2=gint(s1b,s2,1.0D-3,  
     &                    pp,Re,sa,ca,ssig,csig,along1)
                     tec4=tec1+tec2
                     write(16,'(A,I4,1H-,I4,3H),(,I4,1H-,I5,3H),(,
     &I4,1H-,I5,1H))')
     &                    'Electron contents (',nint(h1),nint(h1b),
     &                    nint(h1b),nint(h2),
     &                    nint(h1), nint(h2)
                     write(16,'(16X,3F12.2,A)')
     &                    tec1/1.0D12,tec2/1.0D12,tec4/1.0D12,
     &                    '   x10^15 m^-2'
                  else
                     tec1=gint(s0, s1a,1.0D-3,  pp,Re,sa,ca,ssig,csig,
     &                    along1)
                     tec2=gint(s1a,s1b,1.0D-2,  pp,Re,sa,ca,ssig,csig,
     &                    along1)
                     tec3=gint(s1b,s2, 1.0D-2,  pp,Re,sa,ca,ssig,csig,
     &                    along1)
                     tec4=tec1+tec2+tec3
                     write(16,'(A,2(I4,1H-,I4,3H),(),I4,1H-,I5,3H),(,
     &I4,1H-,I5,1H))')
     &                    'Electron contents (',nint(h0),nint(h1a),
     &                    nint(h1a),nint(h1b),
     &                    nint(h1b),nint(h2),
     &                    nint(h0), nint(h2)
                     write(16,'(16X,4F12.2,A)')
     &                    tec1/1.0D12,tec2/1.0D12,tec3/1.0D12,
     &                    tec4/1.0D12,'   x10^15 m^-2'
                  endif
               endif
            endif
         endif
      endif

      close(16)
      write(6,*)'Output in ',filen1
      end

      real*8 function gint (g1,g2,eps,pp,Re,s1,c1,ssig,csig,along1)
C     special; integrates ELD
      implicit real*8 (a-h,o-z)
      gint1 = 0.0D0
      n = 8
 1    h  =  (g2-g1) / dfloat(n)
      hh = 0.5D0*h
      g  = h*0.5773502691896D0
      y  = g1 + (h-g)*0.5D0
      gint2  =  eld(y,pp,Re,s1,c1,ssig,csig,along1)+
     +     eld(y+g,pp,Re,s1,c1,ssig,csig,along1)
      do m = 1,n-1
         gint2 = gint2 + eld(y+h+g,pp,Re,s1,c1,ssig,csig,along1)+
     +        eld(y+h,pp,Re,s1,c1,ssig,csig,along1)
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
C     special; integrates NeNeQ for vertical profile
      implicit real*8 (a-h,o-z)
      real*8 NeNeQ
      dimension hm(3),aep(3),bb(6)
      common /nepar/hm,aep,bb
      gint1 = 0.0D0
      n = 8
 1    h  =  (g2-g1) / dfloat(n)
      hh = 0.5D0*h
      g  = h*0.5773502691896D0
      y  = g1 + (h-g)*0.5D0
      gint2  =  NeNeQ(y,hm,aep,bb)+
     +     NeNeQ(y+g,hm,aep,bb)
      do m = 1,n-1
         gint2 = gint2 + NeNeQ(y+h+g,hm,aep,bb)+
     +        NeNeQ(y+h,hm,aep,bb)
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
     &     pp,Re,sp,cp,s2,c2,ssig,csig,along1)
      implicit real*8 (a-h,o-z)
      akappa=Re/(Re+400.0D0)
 10   write(6,*)'INPUT: ',
     &     'Ray endpoint 1: latitude (deg N), longitude (deg E),
     &height (km)'
      read(5,*)ph1,alng1,h1
      write(6,*)'INPUT: ',
     &     'Ray endpoint 2: latitude (deg N), longitude (deg E),
     &height (km)'
      read(5,*)ph2,alng2,h2
C     
C     put "vertical" for "near vertical"
      if (abs(ph2-ph1).lt.1.0D-5.and.abs(alng2-alng1).lt.1.0D-5) then
         ph2=ph2
         alng2=alng1
      endif
      r1=Re+h1
      r2=Re+h2
C     
C     provides coordinates (pp,php,alamp) of ray perigee,
C     zenith angle (zeta) for the ray endpoint 1,
C     slant to vertical projection factor (cchi) for the
C     conditions given by akappa (<1)
      call naut(r1,r2,ph1,ph2,alng1,alng2,akappa,
     &     pp,php,alamp,zeta,cchi)
      if (abs(zeta).gt.90.0D0.and.pp.lt.Re) then
         write(6,*) ' ray cuts surface of Earth'
         write(6,*) '      or endpoint 2 lower than endpoint 1.'
         write(6,*) ' Repeat input'
         goto 10
      endif
C     parameters of ray properties for use in other modules (formerly in
C     a common block)
C     (origin for the s coordinate is the ray perigee with
C     lat-long coordinates php and alamp and the radius pp)
C     sp,cp:     sine and cosine of latitude point p (= ray perigee)
C     s2,c2:     sine and cosine of latitude point 2
C     ssig,csig: sine and cosine of ray azimuth
C     (point 2 seen from point 1)
      if (pp.ge.0.1D0) then
         call gcirc(php,ph2,alamp,alng2,sp,cp,s2,c2,ssig,csig,psi)
      endif
C     along1: longitude of point 1 (= ray perigee)
      along1=alamp
      return
      end

      subroutine dat_t_sa(k,iyr,mth,nday,tt,R12,flx)
C     set date and solar activity

C     k=0: input of month only
C     k=1: input of year and month
C     k=2: input of year, month and day of month

C     k=3: input of month and UT
C     k=4: input of year, month and UT
C     k=5: input of year, month, day of month and UT

C     k=6: input of month and LT
C     k=7: input of year, month and LT
C     k=8: input of year, month, day of month and LT

      implicit real*8 (a-h,o-z)
      dimension R12y(12)
C     integer iost
      character*1 yn,fs
      
      open(99,file='R12.dat',status='OLD')
      do while (.true.)
         read(99,*,end=999) iymax
      enddo
  999 continue
      close(99)
      
    1 goto (100,110,120,130,140,150,160,170,180) k+1
 100  iyr=1990
      nday=15
      tt=0.0D0
      write(6,*)'Input: month:'
      read(5,*)mth
      goto 200
 110  nday=15
      tt=0.0D0
      write(6,*)'Input: year, month:'
      read(5,*)iyr,mth
      goto 200
 120  tt=0.0D0
      write(6,*)'Input: year, month, day of month:'
      read(5,*)iyr,mth,nday
      goto 200
 130  iyr=1990
      nday=15
      write(6,*)'Input: month, UT:'
      read(5,*)mth,tt
      goto 200
 140  nday=15
      write(6,*)'Input: year, month, UT:'
      read(5,*)iyr,mth,tt
      goto 200
 150  write(6,*)'Input: year, month, day of month, UT:'
      read(5,*)iyr,mth,nday,tt
      goto 200
 160  write(6,*)'Input: month, LT:'
      read(5,*)mth,tt
      goto 200
 170  write(6,*)'Input: year, month, LT:'
      read(5,*)iyr,mth,tt
      goto 200
 180  write(6,*)'Input: year, month, day of month, LT:'
      read(5,*)iyr,mth,nday,tt
      goto 200
  200 continue
      if (iyr.gt.100.and.(iyr.lt.1931.or.iyr.gt.2049)
     &     .or.iyr.lt.0) then
         write(6,*)
     &        'error in year (valid: 1931-2049 or 0-49 for 2000-2049'
         write(6,*)'   or 50-99 for 1950-1999)'
         write(6,*)'Repeat'
         goto 1
      endif
      if (mth.lt.1.or.mth.gt.12.or.tt.lt.0.0D0.or.tt.gt.24.0D0) then
         write(6,*)
     &        'input of month or UT not valid (valid: 1-12 and 0-24)'
         write(6,*)'Repeat'
         goto 1
      endif
      if (iyr.lt.50) iyr=iyr+2000
      if (iyr.lt.1900) iyr=iyr+1900
      if (k.ne.0.and.k.ne.3.and.iyr.ge.1931.and.iyr.le.iymax) then
         write(6,*)'User input R12/F10.7 for this year and month? (y/n)'
         read(5,'(A)')yn
      else
         yn='y'
      endif
      if (yn.eq.'n'.or.yn.eq.'N') then
         open(99,file='R12.dat',status='OLD')
 211     read(99,*)j,R12y
         if (j.lt.iyr) goto 211
         close(99)
         R12=R12y(mth)
         if (R12.gt.150.0D0)then
            write(6,'(A/2A/A)')'***WARNING!Sunspot number exceeds',
     &           ' 150*** R12=150 will be used',
     &           ' (ITU-R P.1239 recommendation).'
            R12=150.0D0
         endif
         flx=63.7D0+(0.728D0+8.9D-4*R12)*R12
      else
 214     write(6,*)
     &        'INPUT: solar activity type:',
     &        ' sunspot number (S) or 10.7 cm radio flux (F)?'
         read(5,'(A)')fs
         if (fs.eq.'F'.or.fs.eq.'f') then
 212        write(6,*)'INPUT: radio flux (>=63 units)'
            read(5,*)flx
            if (flx.lt.63.0D0) then
               write(6,*)'***WARNING!Solar flux F<63 F.U. might lead'
               write(6,*)'to undefined electron density values'
               goto 212
            endif
            if (flx.gt.193.0D0)then
               write(6,'(A/2A)')'***WARNING!Solar flux F exceeds ',
     &              ' 193 F.U.*** F= 193 F.U. will be used',
     &              ' (ITU-R P.1239 recommendation).'
               flx=193.0D0
            endif
            R12=sqrt(167273.0D0+(flx-63.7D0)*1123.6D0)-408.99D0
         else
            if (fs.eq.'S'.or.fs.eq.'s') then
 213           write(6,*)'INPUT: sunspot number (R12)'
               read(5,*)R12
               if (R12.lt.0D0) then
                  write(6,*)'***WARNING: R12<0 might lead'
                  write(6,*)'to undefined electron density values'
                  goto 213
               endif
               if (R12.gt.150.0D0)then
                  write(6,'(A/2A)')'***WARNING!Sunspot number S ',
     &                 ' exceeds 150*** R12=150 will be used',
     &                 ' (ITU-R P.1239 recommendation).'
                  R12=150.0D0
               endif
            else
               goto 214
            endif
            flx=63.7D0+(0.728D0+8.9D-4*R12)*R12
         endif
      endif
      return
      end
      
      subroutine gcirc(alat1,alat2,along1,along2,s1,c1,s2,c2,ssig,
     &     csig,psi)
C     
C     calculates great circle path properties
C     input parameters are the endpoints coordinates
C     (alat1,along1), (alat2,along2) [deg.]
C     output parameters are sine and cosine of the endpoint latitudes
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
C     gives electron density as a function of the coordinate s
C     (distance from the perigee of the ray)
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
