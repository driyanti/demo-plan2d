!*******************************************************************
! name:   obsmodel
!
! descr:  
!
! NOTES:
!
! HIST:   Version Date     Author       Comment
!         1.0     09/09/17 rgvb          Initial version
!         2.0     02/03/20 christina     Modified
!*******************************************************************

!**********************************************************************
!    declare variables
!**********************************************************************
program plan1d
      
implicit none
include "mathConstants.h"

integer  :: kk, nf,nmult, nf1,nf2
integer  :: nlay,i,j,k,l,ntrac,nsamp,nwav,iom
integer  :: nx, direct
real     :: df
real     :: dw,dk,dkx,dl
real     :: f1,f2   
real     :: dx,dt,x3r,x3s,eps

real, dimension(nmax)       :: ktap, kx2
real, dimension(nmax)       :: wave, etap
real, dimension(nmax)       :: velo, rho, depth
real, dimension(nmax, nmax) :: pupdat, pdat, vdodat

complex ::    cdirct,cwsurf,csghos,crghos,cprim,cmult1,crho,crho_end,ctt
complex ::    cmultt,cgtot,cmults,crlay,coef,ck,ces

complex, dimension(10, nmax)   :: ck2l
complex, dimension(nmax)       :: csgamma, cpup, cpdo, cvdo
complex, dimension(nmax)       :: cp, cpr, cdirect, cwsurfr, &
                                  cr, cwave, ct, ct1, ct2, cwave2
complex, dimension(nmax, nmax) :: cpupdat, cvdodat, cpdat 

!*********************************************************************C
!    define constants       
!*********************************************************************C
dt    = 0.004
nf    = 2048 
nsamp = 1024
nx    = 1024
ntrac = 501
x3s   = 5.
x3r   = 25.
dx    = 10. 
df    = 1./(nf*dt)
f1    = 2.
f2    = 80.
nf1   = int (f1/df)
nf2   = int (f2/df)
dw    = 2.*pi*df
dl    = 1./(nx*dx)
dkx   = 2.*pi*dl
eps   = 0.5  
nmult = 4
direct = 2

!**********************************************************************
!    read the geometry from the ascii file                     
!**********************************************************************

open (10,file='/home/dwi/Dwi/geo4',status='unknown')
read(10,*)nlay
do i=1,nlay
   read(10,*)depth(i),velo(i), rho(i)
enddo
close(10)

!*********************************************************************C
!    calculate (omega/v)**2 for all omega
!*********************************************************************C

do l = 1,nlay   
   dk = dw / velo(l)
   ck = cmplx (eps / velo(l), -dk)
   do iom = 1, nf/2+1   
      ck          = ck + ci * dk
      ck2l(l,iom) = ck ** 2
   enddo  
enddo

write(*,*)'ready with zero-part ...', nmax, ck2l(1,100)

!*********************************************************************C
!    calculate dk**2 for all kx
!*********************************************************************C

do j = 1,nx/2+1
   kx2(j) = ((j-1)*dkx)**2
enddo  

!*********************************************************************C
!    calculate the time-taper    
!*********************************************************************C

call timetap(etap,nsamp,dt,eps)

!*********************************************************************C
!    calculate the spacial tapers    
!*********************************************************************C

call gen_taper(ktap,nx,1,1,nx/2-nx/8,nx/2,2.0)

!*********************************************************************C
!    calculate the wavelet spectrum
!*********************************************************************C
open(15,file='/home/dwi/Dwi/wave.save',status='unknown')
read(15,*)nwav

do i=1,nwav
   read(15,*)wave(i) 
enddo      
close(15)

cwave(1:nf*2) = cmplx(0.0,0.0)
!call vclr(cwave, nf*2)        

do i=1,nwav  
   cwave(i)=cmplx (wave(i), 0.)
enddo 

write(*,*)'ready with cwave ...ok', cwave(100)

call pfft(cwave,nf,-dt)

!    ... use cwave2 below to get an artifact-free vdown field ...
!    ... use cwave if you want the best pdemult result ...

do iom=1,nf/2
   !cwave2(iom)  = cmplx (eps, (iom-1) * dw) * rho(1) * cwave(iom)
   cwave2(iom) = cmplx (1, 0.) * cwave(iom)
enddo

write(*,*)'ready with cwave2 ...', rho(1), cwave2(100)
!*********************************************************************C
!    start of frequency loop 
!*********************************************************************C
write(*,*)'ready with frequency loop ...', nmax

!... scale array to zero ...        
cp(1:nx*2)      = c0
cpup(1:nx*2)    = c0
cpdo(1:nx*2)    = c0
cvdo(1:nx*2)    = c0
cdirect(1:nx*2) = c0
cwsurfr(1:nx*2) = c0
cpr(1:nx*2)     = c0

do iom=nf1,nf2              

    !... calculate constants ...
    ces = cmplx(eps,(iom-1)*dw)

    !nkmax=int((iom-1)*dw / (1500.*dkx)) 
    !call gen_taper(ktap,nx,1,1,nkmax-nkmax/8,nkmax,2.0)

    !... loop for all wavenumbers ...
    do j=1,nx/2+1

        !... calculate the gamma's and the reflection coefficients ...
        csgamma (nlay) = csqrt (ck2l(nlay,iom) + kx2 (j))
        crho_end = csgamma (nlay) / rho (nlay)
        cr (nlay) = 0.
        do k = nlay - 1, 1, - 1
            csgamma(k) = csqrt (ck2l(k,iom) + kx2 (j))
            crho = csgamma (k) / rho (k)
            ctt = (crho - crho_end) / (crho + crho_end)
            crlay = cr (k + 1) * cexp (-2. * csgamma (k + 1) &
                  * (depth (k + 1) - depth (k)))
            cr (k) = (ctt + crlay) / (1. + ctt * crlay)
            crho_end = crho
        enddo

        coef = cr (1)

        !... calculate the direct wave ...
        if (direct==1) then
            cdirct = cexp (-csgamma(1) * abs (x3s - x3r))
        else
            cdirct = 0.0
        endif

        !... calculate the water surface reflection ...
        cwsurf = -1. * cexp (-csgamma(1) * (x3s + x3r))

        !... calculate the source ghost delay factor ...
        csghos = 1.  - cexp ( -2. * csgamma(1) * x3s)

        !... calculate the receiver ghost delay factor ...
        crghos = 1. - cexp ( -2. * csgamma(1) * x3r)

        !... calculate the primary reflection ...
        cprim  = cexp (-csgamma(1) * (2.*depth(1) - x3s - x3r))

        !... calculate the multiples ...
        cmult1 = -1. * coef * cexp ( -2. * csgamma(1) * depth(1))
        cmults = cmult1 ** (nmult + 1)
        cmultt = (1. - cmults) / (1. - cmult1)

        !... calculate the total green's function for the pressure field ...
        cgtot = cdirct + cwsurf + coef * csghos * crghos * &
                cprim * cmultt

        !... compute the final pressure wavefield ...
        cp(j)      = cwave2(iom) * cgtot  / (2. * csgamma(1))
        cdirect(j) = cwave2(iom) * cdirct / (2. * csgamma(1))
        cwsurfr(j) = cwave2(iom) * cwsurf / (2. * csgamma(1)) 
        cpr(j)     = cwave2(iom) * cprim / (2. * csgamma(1))

        !... compute up and down going wavefields (at x3=x3r) ...
        cpup (j) = (cp(j)-(cdirect(j)+cwsurfr(j))) / crghos
        cpdo (j) = cp(j) - cpup(j)

        !... compute downgoing vz using 10.31 ...
        cvdo (j) =  cpdo (j) * csgamma(1) / (ces * rho(1))

        !... apply tapers ... 
        cp (j)   = cp (j)   * ktap(j)
        cpup (j) = cpup (j) * ktap(j)
        cpdo (j) = cpdo (j) * ktap(j)
        cvdo (j) = cvdo (j) * ktap(j)
        cpr  (j) = cpr (j) * ktap(j)
    
        !... make symmetrical for the negative kx values ...
        cp (nx + 2 - j)   = cp (j)
        cpup (nx + 2 - j) = cpup (j)
        cpdo (nx + 2 - j) = cpdo (j)
        cvdo (nx + 2 - j) = cvdo (j)
        cpr (nx + 2 - j)  = cpr (j)
    enddo

    !... transform back to space ...cpr, cpup, cvdo 

    call pfft(cpr,nx,dl)    
    call pfft(cpup,nx,dl)
    call pfft(cvdo,nx,dl)

    do j = 1,ntrac/2
        cpdat (j,iom) = cp (nx-ntrac/2+j)
        cpupdat (j,iom) = cpup (nx-ntrac/2+j)
        cvdodat (j,iom) = cvdo (nx-ntrac/2+j)
    enddo

    do j = ntrac/2+1,ntrac
        cpdat (j,iom) = cpr (j-ntrac/2)
        cpupdat (j,iom) = cpup (j-ntrac/2)
        cvdodat (j,iom) = cvdo (j-ntrac/2)
    enddo

!    ... end of frequency loop ...
enddo

write(*,*)'ready with the model...', cpdat(100,100)

!*********************************************************************C
!    go back to time and write files to disk 
!*********************************************************************C

do kk=1,ntrac
    !... scale array to zero ...
    call vclr(ct,nf*2)        
    call vclr(ct1,nf*2)
    call vclr(ct2,nf*2)

    !... set negative frequencies ...
    do iom = 1, nf/2+1  
        ct(iom)  = cpdat(kk,iom)
        ct1(iom) = cvdodat(kk,iom)
        ct2(iom) = cpupdat(kk,iom)
    enddo

    !... transform back to time ...
    call pfft(ct,nf,df)
    call pfft(ct1,nf,df)
    call pfft(ct2,nf,df)

    !.... apply taper ...
    do i = 1,nf/2 
        pdat(kk,i)   = 2. * real(ct(i))
        vdodat(kk,i) = 2. * real(ct1(i))
        pupdat(kk,i) = 2. * real(ct2(i))
    enddo

!tel = tel+1
!CALL writedat(lunpsct,pdat(kk,i),nsamp,tel)

enddo
!close(lunpsct)

write(*,*)'ready with the model...dpr', pdat(100,100)
OPEN(unit=75, &
     file='/home/dwi/Output/pdat.b',&
     access='direct', form='unformatted',recl=4*nsamp*ntrac)

write(75,rec=1)((pdat(i,kk),kk=1,nsamp),i=1,ntrac)
close(75)

write(*,*)'ready with the model...dpr'
OPEN(unit=70,file='/home/dwi/Output/pup.b', &
     access='direct', form='unformatted',recl=4*nsamp*ntrac)
write(70,rec=1)((pupdat(i,kk),kk=1,nsamp),i=1,ntrac)  

OPEN(unit=80,file='/home/dwi/Output/vdo.b', &
     access='direct', form='unformatted',recl=4*nsamp*ntrac)
write(80,rec=1)((vdodat(i,kk),kk=1,nsamp),i=1,ntrac)

close(70)
close(80)

!tel = tel+1
!CALL writedat(lunpsct,pdat(kk,i),nsamp,tel)
write(*,*)'ready with the model of pup and vdo-parts ...'

!close(lunpsct) 
!*********************************************************************c
!    end of program    
!*********************************************************************c
return 
end
!**********************************************************************
!    close all files

!**********************************************************************
!    some general utilities
!**********************************************************************

      subroutine timetap(t,n,dt,eps)

      implicit none

      integer i,n
      real    eps,dt,t(n)

!    .... calculate the time taper ....
      do i = 1,n
         t(i) = exp ((i-1)*dt*eps)
      enddo

!    ... end of subroutine ...
      return
      end
!**********************************************************************
      subroutine vclr(ct,n)
      
      integer i, n
      complex ct(n)
      
      do i=1,n
       ct(i) = cmplx(0.0,0.0)
      end do
      
      return
      end

SUBROUTINE OPENDA(LUN,FILNAM,NBYTES)
 
!C     ... DECLARATION OF PARAMETERS ...

	  IMPLICIT NONE
 
      INTEGER LUN,IOS,NBYTES
      CHARACTER*(*) FILNAM

!C     ... BEGIN OF SUBROUTINE ...

      OPEN (UNIT   = LUN, FILE = FILNAM, FORM   ='UNFORMATTED',&
            ACCESS ='DIRECT', RECL   = NBYTES, ERR = 999, IOSTAT = IOS)
	               
!C     ... REWIND THE FILE ...
      REWIND(LUN)

      RETURN

!C     ... ERROR HANDLING ...
999	  WRITE (*,*) 'ERROR OCCURED IN OPENING FILE ON UNIT NUMBER',LUN
	  WRITE (*,*) 'THE NAME OF THE FILE IS',FILNAM
	  WRITE (*,*) 'THE ERROR CODE IS',IOS

!C     ... END OF SUBROUTINE ...
	  RETURN
      END
!C----------------------------------------------------------------------
!**********************************************************************

      SUBROUTINE WRITEDA (LUN,DAT,NSAMP,RECNR)

!HEADER**********************************************************
!TITLE        : WRITEDA
!DESCRIPTION  : THIS SUBROUTINE WRITES NSAMP IN A DIRECT ACCESS FILE.
!CK  KEYWORDS : WRITEDA
!C*
!CL  LANGUAGE     : FORTRAN77

!CALLING SEQ. : CALL WRITEDA(LUN,DAT,NSAMP,RECNR)

!INPUTPAR.    : LUN    - LOGICAL UNIT NUMBER OF FILE
!CI               : DAT    - ARRAY OF DATA TO BE FILLED
!CI               : NSAMP  - NUMBER OF REAL SAMPLES TO BE READ
!CI               : RECNR  - RECORDNUMBER TO BE READ
!CO  OUTPUTPAR.   : DAT    - ARRAY WITH DATA
!ENTRY POINTS     : NONE
!CB  COMMON BLOCK : NONE
!***********************************************HEADERC
 
!C     ... DECLARATION OF VARIABLES ...

      IMPLICIT NONE 

	  INTEGER RECNR,LUN,NSAMP   
	  REAL DAT(NSAMP)
 
!C     ... READ THE DATA ...
	  WRITE (LUN,REC=RECNR,ERR=10) DAT 
 
	  RETURN

!C     ... ERROR HANDLING ...
10	  WRITE(*,*)'ERROR OCCURRED WHILE WRITING RECORD IN WRITEDA &
                FROM UNIT',LUN
!C     ... END OF SUBROUTINE ...
	  RETURN
      END


!**********************************************************************

      subroutine pfft(cx,n,delta)

!    ... implicit none ...

      integer j,i,N,m,istep,l
      complex cx(N),cw,ctemp
      real signi,sc,delta,arg

!    ... start of subroutine ...
      j=1
      if (delta.lt.0.) then
          signi     = -1.
      else
          signi     = 1.
      endif
      sc  = abs(delta)
      do 630 i=1,N
      if(i.gt.j) goto 610
      ctemp=cx(j)*sc
      cx(j)=cx(i)*sc
      cx(i)=ctemp
610   m=N/2
620   if(j.le.m) goto 630
      j=j-m
      m=m/2
      if(m.ge.1) goto 620
630   j=j+m
      l=1
640   istep=2*l
      do 650 m=1,l
      arg=3.141592653*signi*(m-1)/l
      cw=cmplx(cos(arg),sin(arg))
      do 650 i=m,N,istep
      ctemp=cw*cx(i+l)
      cx(i+l)=cx(i)-ctemp
650   cx(i)=cx(i)+ctemp
      l=istep
      if(l.lt.N) goto 640

!    ... end of subroutine ...
      return
      end

!    create cos**power taper from itap1 to itap2 and itap3 to itap4
      subroutine gen_taper(taper,lentap,itap1,itap2,itap3,itap4,power)
 
      implicit none

!**********************************************************************
!    subroutine variables
!**********************************************************************

      integer itap1,itap2,itap3,itap4,lentap,ix
      real power
      real taper(*)

!**********************************************************************
!    local variables
!**********************************************************************

      real tapfac,pi

!**********************************************************************
!    define pi
!**********************************************************************

      pi=3.141592654

!**********************************************************************
!    check input variables
!**********************************************************************

      if (lentap.le.0) then
         write(*,*)'length taper < 1'
      endif
      if ( itap2.lt.itap1 .or. itap3.lt.itap2 .or. itap4.lt.itap3 ) then
         write(*,*)'gen_taper: taper values should increase'
         write(*,*)'wrong taper values!'
      endif

!**********************************************************************
!    create taper, start with zero part until itap1
!**********************************************************************

      do ix=1,itap1-1
         taper(ix)=0.
      enddo

!**********************************************************************
!    cosine edge from itap1 to itap2
!**********************************************************************

      tapfac=0.5*pi/float(max(1,itap2-itap1+1))
      do ix=itap1,itap2
         taper(ix)=abs(cos(float(ix-itap2)*tapfac))**power
      enddo

!**********************************************************************
!    flat part with value 1 from itap2 to itap3
!**********************************************************************

      do ix=itap2+1,itap3-1
         taper(ix)=1.
      enddo

!**********************************************************************
!    cosine edge from itap3 to itap4
!**********************************************************************

      tapfac=0.5*pi/float(max(1,itap4-itap3+1))
      do ix=itap3,itap4
         taper(ix)=abs(cos(float(ix-itap3)*tapfac))**power
      enddo

!**********************************************************************
!    zero beyond itap4
!**********************************************************************

      do ix=itap4+1,lentap
         taper(ix)=0.
      enddo

!**********************************************************************
!    end of subroutine GEN_TAPER
!**********************************************************************

      return
      end

