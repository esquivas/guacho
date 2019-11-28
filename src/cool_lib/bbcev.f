#ifdef COOLINGBBC
!=======================================================================
!   Cooling routine using Benjamin, Bensonm, & Cox (2001) ApJL
!   non-equilibrium ionization routines.
!=======================================================================
      subroutine bbcev(deltat,denH,TTin,yi,y)
c========================================================================
c      INPUT:  Timestep (seconds) deltat
c              Hydrogen density (cm-3)           
c              Temperature (K) TTin
c              Input Ionization vector (see README.bbc for details)
c     --------------------
c     NOTE: this vector includes(!) abundances of ions.
c     For example,y(He I)=y(3) ranges from 0 to 0.1!
c     y(1)=f(H I)
c     y(2)=f(H II)
c     y(3)=f(He I)*abund(He)
c     y(4)=f(He II)*abund(He)
c     y(5)=f(He III)*abund(He)
c     y(NTOT-1)=zbar (mean charge of trace ions)
c     y(NTOT)=electron fraction
c     --------------------
c     OUTPUT: Time evolved ionization vector
c     Written by Bob Benjamin 20 Jan 2001
c------------------------------------------------------------------
#include "bbc.h"
      dimension yi(*),y(*)
      dimension crec(NIONS),cion(NIONS)

      parameter(NWORK=30)
      dimension A(30),B(30),C(30),f(30)

c     ---------------------------
c     Copy the initial ionization
c     fraction array into the 
c     working array
c     ---------------------------
      do 5 i=1,NTOT
 5       y(i)=yi(i)  


c     --------------------------
c     Zero out rate coefficient
c     arrays c     -------------------------
      do 55 i=1,NIONS
         crec(i)=0.
         cion(i)=0.
 55   continue

c     -----------------------------------
c     Check to make sure temperature lies
c     in array bounds
c     -----------------------------------
      Tlog=alog10(TTin)
      if ( (Tlog.lt.T_(1)).or.(Tlog.gt.T_(NTEMP)) ) then 
        print *, 'T is out of range...', TTin
        !TTin=1.e4
        !Tlog=4.0
        stop
      endif

c     ----------------------------------
c     iT is index such that Tlog lies
c     between T_(iT) and T_(iT+1)
c     wl and wu are coefficients for
c     linear interpolation
c     ----------------------------------
      dTlog=(T_(NTEMP)-T_(1))/(NTEMP-1)
      iT=int( (Tlog-T_(1))/dTlog )+1
      wl=( 10**(T_(iT+1))-TTin )/(10**(T_(iT+1))-10**(T_(iT)))
      wu= 1 - wl

c     ----------------------------------
c     Loop backwards through elements to
c     calculate ionization balance
c     ----------------------------------
 1000 continue

      dene=denH*y(NTOT)
c     ----------------------------------
c     Fill up array of rate tables, from
c     top down. 
c     ----------------------------------
      
      ntop=NTOT-(NOTHER+1)
      do 200 ii=ntop,1,-1

c        radiative recombination
         crec(ii)=dene*(wl*rrec(ii,iT)+wu*rrec(ii,iT+1))

c        collisional ionization
         cion(ii)=dene*(wl*rcol(ii,iT)+wu*rcol(ii,iT+1))

 200  continue

      do 100 iel=1,NEL
         ix=indexel(iel)
         ilen=noel(iel)+1

         big=0.

         do 101 i=1,ilen

            tmp=cion(ix+i-1)+crec(ix+i-1)
            if ((yi(ix+i-1)/abund(ix+i-1)).gt.1e-4) big=max(big,tmp)
 101     continue

         jsub=deltat*big
         jsub=max(1,min(jsub,100))
         if(iel.eq.1) then
            jsub=1
            dt=deltat
         else
            dt=deltat/jsub
         end if


         do 199 j=1,jsub
           if(iel.eq.1) then
              f(1)=y(1)/abund(1)
              f(1)=f(1)*exp(-(crec(2)+cion(1))*dt)+crec(2)/
     >              (cion(1)+crec(2))*(1-exp(-(crec(2)+cion(1))*dt))
              f(2)=1.-f(1)
              y(1)=f(1)*abund(1)
              y(2)=f(2)*abund(2)
           else
              do 113 i=1,ilen
                f(i)=0.
                A(i)=0.
                B(i)=0.
                C(i)=0.
 113          continue
              do 102 i=1,ilen
                 f(i)  = y(ix+i-1)/abund(ix+i-1)
                 B(i)  = 1+(cion(ix+i-1)+crec(ix+i-1))*dt
 102          continue
              do 112 i=1,ilen-1
                 A(i+1)= -cion(ix+i-1)*dt
                 C(i)  = -crec(ix+i)*dt
 112          continue
              call tridag_inv(A,B,C,f,f,ilen)
              do 103 i=1,ilen
 103             y(ix+i-1)=f(i)*abund(ix+i-1)
           end if
            
 199    continue
 100  continue

c     ------------------------------------
c     Evolve zbar...
c     ------------------------------------
      zbarnew=zbarneqr(dene,TTin,y(NTOT-1),deltat)
      y(NTOT-1)=zbarnew
c     -------------------------------------
c     Recalculate electron density. 
c     -------------------------------------
      y(NTOT)=yelect(y)+(1.0e-9+0.00144313*y(NTOT-1))

      return
      end
c=======================================================================BBCEV


       function bbcD(dene,TTin,ZZin)
c========================================================================BBCD
c      INPUT: electron density (cm-3)
c             Temperature (K)
c             Mean charge state (Zbar) [ranges from 0.1833 to 8.377]
c      OUTPUT: Rate of change of Zbar D(Z,T) in s-1. 
c      Written by Bob Benjamin 20 Jan 2001
c       include "bbc.h"
#include "bbc.h"
       tolerance=0.05

       Tin=alog10(TTin)
       Zin=ZZin       


c      If input parameters lie off the precomputed
c      grid, use ionization/recombination rate at
c      the edge of the grid...

       if (Tin.lt.Tneq(1)) then 
         Tin=Tneq(1)*1.001
       else if (Tin.gt.Tneq(NTneq)) then 
         Tin=Tneq(NTneq)*0.999
       endif

       if (Zin.lt.Zneq(1)) then 
         Zin=Zneq(1)*1.001
       else if (Zin.gt.Zneq(NZneq)) then 
         Zin=Zneq(NZneq)*0.999
       endif


       call hunter(Tneq,NTneq,Tin,ilo)
       call hunter(Zneq,NZneq,Zin,jlo)
       si=(Tin-Tneq(ilo))/(Tneq(ilo+1)-Tneq(ilo))
       tj=(Zin-Zneq(jlo))/(Zneq(jlo+1)-Zneq(jlo))
       

       if ((ilo.eq.0).or.(ilo.eq.NTneq)) then
         print *, 'bbcD: temperature out of bounds', Tin
c         stop       ! Aug012003-jc
         bbcD=0.     ! Aug012003-jc
         go to 111   ! Aug012003-jc
       endif  
       if ((jlo.eq.0).or.(jlo.eq.NZneq)) then 
          
         print *, 'bbcD: ionization level out of bounds', Zin
c         stop
         bbcD=0.
         go to 111
       endif

       ZTeq=(1-si)*Zeq(ilo)+si*Zeq(ilo+1)
       

       dZ=Zin-ZTeq
c      Determine if we are essentially converged..
       if (abs(dZ).LT. 0.002) then
        bbcD=0.
        return
       endif

c      Determine whether ionizing or recombining
       if (dZ .LT. 0) then 
        iion=1
       else
        iion=-1
       endif
      
c      Determine whether we should use a fine grid
       if (abs(dZ).LT. 0.9) then
        ifine=1
       else 
        ifine=0
       endif     

c     Calculate the time derivative of recom rate...
      if (ifine.eq.0) then 
         Din=(1-si)*(1-tj) *alog10(abs(D(ilo  , jlo  )))+
     >          si  *(1-tj)*alog10(abs(D(ilo+1, jlo  )))+
     >        (1-si)*  tj  *alog10(abs(D(ilo  , jlo+1)))+        
     >          si  *  tj  *alog10(abs(D(ilo+1, jlo+1)))
      else 
         Zd=abs(dZ)
         call hunter(Zf,NZneqf,Zd,jlo)
         tj=(Zd-Zf(jlo))/(Zf(jlo+1)-Zf(jlo))       
c        Fudge 1
         if (iion.eq.1) then 
           if (ilo.eq.1) then 
             si=0.
           else
             ilo=ilo-1
           endif
           Din=(1-si)*(1-tj) *alog10(abs(Di(ilo  , jlo  )))+
     >            si  *(1-tj)*alog10(abs(Di(ilo+1, jlo  )))+
     >          (1-si)*  tj  *alog10(abs(Di(ilo  , jlo+1)))+        
     >            si  *  tj  *alog10(abs(Di(ilo+1, jlo+1)))

         else 
c          Fudge 2
           if (ilo.eq.NTneq-1) then 
             ilo=ilo-1
             si=1.0
           endif
           Din=(1-si)*(1-tj) *alog10(abs(Dr(ilo  , jlo  )))+
     >            si  *(1-tj)*alog10(abs(Dr(ilo+1, jlo  )))+
     >          (1-si)*  tj  *alog10(abs(Dr(ilo  , jlo+1)))+        
     >            si  *  tj  *alog10(abs(Dr(ilo+1, jlo+1)))
         endif
       endif
       Dinlog=Din
       dzdf=iion*10**(Din)

c      Note that the fudge1/fudge2 involve setting the recom/ion
c      rate to the value at either 4.1 or 7.9 if the temperature
c      lies in the domain corner...(e.g., for region between 10^4 
c      and 10^4.1 there is no left side for ionizing gas...)
c      In addition, remember that Di[1,NTneq-1] ranges from 4.1 to 8.0
c      Dr[1,NTneq-1] ranges from 4.0 to 7.9

       bbcD=dene*dzdf

 111   continue
       return
       end
c========================================================================BBCD

c=======================================================================
      subroutine  tridag_inv(a,b,c,f1,f2,ilen)
c=======================================================================
      parameter (nmax=100)
      dimension gam(nmax),a(ilen),b(ilen),c(ilen)
      dimension f1(ilen),f2(ilen)

      b1=b(1)
      f2(1)=f1(1)/b1
      do i=2,ilen
        gam(i)=c(i-1)/b1
        b1=b(i)-a(i)*gam(i)
        if (b1.EQ.0.) pause
        f2(i)=(f1(i)-a(i)*f2(i-1))/b1
      enddo

      do  j=ilen-1,1,-1
        f2(j)=f2(j)-gam(j+1)*f2(j+1)
      enddo

      return
      end
c=======================================================================

c=======================================================================
      function zbarneqr(dene,Tk,zinit,dt)
c=======================================================================
c      include "bbc.h"
#include "bbc.h"

c     calculation ionization/recomb time (in scaled units)
      dzdt=bbcD(dene,TK,zinit)

c     exit if no zbar evolution      
      if (dzdt.eq.0) then 
        zbarneqr=zinit
        return
      endif

c     calculate which timestep is shorter, 
c     the input timestep, or 1% of the 
c     ionization time step
      tz=0.01*zinit/abs(dzdt)
      dtmin=min(dt,tz)
      NTZ=4*int(dt/dtmin)
      NTZ=max(5,NTZ)
      dtz=dt/NTZ
 
      z0=zinit
      do iz=1,NTZ 
        z0=z0+bbcD(dene,Tk,z0)*dtz
      enddo
      zbarneqr=z0

      return
      end
c=======================================================================


c=======================================================================BBCAB
       function bbcab(n_el)
c=======================================================================
c      INPUT: n_el=element atomic number
c      OUTPUT: solar abundance (linear) from Grevesse \& Anders 
c      Written by Bob Benjamin 20 Jan 2001
c------------------------------------------------------------------------
       if (n_el.eq.2) then 
           bbcab=0.0977
       else if (n_el.eq.6) then
           bbcab=3.55e-4
       else if (n_el.eq.7) then
           bbcab=9.93e-5
       else if (n_el.eq.8) then
           bbcab=7.41e-4
       else if (n_el.eq.10) then
           bbcab=1.20e-4
       else if (n_el.eq.12) then
           bbcab=3.80e-5
       else if (n_el.eq.14) then
           bbcab=3.55e-5
       else if (n_el.eq.16) then
           bbcab=1.55e-5
       else if (n_el.eq.18) then
           bbcab=3.16e-6
       else if (n_el.eq.20) then
           bbcab=2.29e-6
       else if (n_el.eq.26) then
           bbcab=3.16e-5
       else if (n_el.eq.28) then
           bbcab=1.78e-6
       else
           print *, 'Invalid element number'
           stop
       endif

       return
       end
c========================================================================BBCAB

c     Pasted Aug 12 2002 - AE
c========================================================================BBCCOOL
       function bbccool(zMet,denH,TTin,y)
c========================================================================
c      INPUT:  Metallicity (1=solar/linear units) Zmet
c              Hydrogen particle density (cm-3) denH
c              Temperature (K) TTin
c              Ionization vector (see README.bbc) for details...
c      OUTPUT: Cooling rate for both H/He and trace elements...
c              in units of ergs s-1 cm-3 (n_e*n_H has already been 
c              multiplied in...)
c      Written by Bob Benjamin 20 Jan 2001
      include "bbc.h"
      dimension y(*)

c     ----------------------------
c     calculate electron fraction
c     ----------------------------
      dene=denH*y(NTOT)
      ZZin=y(NTOT-1)

c     ------------------------------------------------------------
c     Calculate hydrogen/helium cooling.....
c     ------------------------------------------------------------

c     -----------------------------------
c     Check to make sure temperature lies
c     in array bounds
c     -----------------------------------
      Tlog=alog10(TTin)
      if ( (Tlog.lt.T_(1)).or.(Tlog.gt.T_(NTEMP)) ) then
       print *, 'Tlog out of range in bbccool.f'
       stop
      endif

c     ----------------------------------
c     iT is index such that Tlog lies
c     between T_(iT) and T_(iT+1)
c     wl and wu are coefficients for
c     linear interpolation
c     ----------------------------------
      dTlog=(T_(NTEMP)-T_(1))/(NTEMP-1)
      iT=int( (Tlog-T_(1))/dTlog )+1
      wl=( 10**(T_(iT+1))-TTin )/(10**(T_(iT+1))-10**(T_(iT)))
      wu= 1 - wl

c     --------------------------------------
c     Calculate the cooling 
c     in units of ergs s-1 
c     The factor of 1e-23 is since the
c     cooling tables are expected to be
c     multiplied by 1e23
c     --------------------------------------
      cool=0.
      do 300 ii=1,NIONS
         cool=cool+y(ii)*(wl*rcool(ii,iT)+wu*rcool(ii,iT+1))
 300  continue

      bbccool=cool*dene*denH*1e-23

c     ------------------------------------------------------------
c     Now calculate the nonequilbrium cooling...
c     ------------------------------------------------------------

       Tin=alog10(TTin)
       Zin=ZZin

c      If input parameters lie off the precomputed
c      grid, use ionization/recombination rate at
c      the edge of the grid...

       if (Tin.lt.Tneq(1)) then 
         Tin=Tneq(1)*1.001
       else if (Tin.gt.Tneq(NTneq)) then 
         Tin=Tneq(NTneq)*0.999
       endif

       if (Zin.lt.Zneq(1)) then 
         Zin=Zneq(1)*1.001
       else if (Zin.gt.Zneq(NZneq)) then 
         Zin=Zneq(NZneq)*0.999
       endif


       call hunter(Tneq,NTneq,Tin,ilo)
       call hunter(Zneq,NZneq,Zin,jlo)
       si=(Tin-Tneq(ilo))/(Tneq(ilo+1)-Tneq(ilo))
       tj=(Zin-Zneq(jlo))/(Zneq(jlo+1)-Zneq(jlo))


       if ((ilo.eq.0).or.(ilo.eq.NTneq)) then
         print *, 'bbccool: temperature out of bounds', Tin
       endif
       if ((jlo.eq.0).or.(jlo.eq.NZneq)) then 
         print *, 'bbccool: ionization level out of bounds', Zin
       endif

c      Check to make sure interpolation constants are
c      in bounds    

c       ZTeq=(1-si)*Zeq(ilo)+si*Zeq(ilo+1)
c       rLTeq=(1-si)*rLeq(ilo)+si*rLeq(ilo+1)


c       dZ=Zin-ZTeq
c      Determine if we are essentially converged..
c       if (abs(dZ).LT. 0.02) then
c        cooleq=rLin
c       endif

c     Interpolate the cooling
         rLin=(1-si)*(1-tj)*alog10(rL(ilo  , jlo  ))+
     >          si  *(1-tj)*alog10(rL(ilo+1, jlo  ))+
     >        (1-si)*  tj  *alog10(rL(ilo  , jlo+1))+        
     >          si  *  tj  *alog10(rL(ilo+1, jlo+1))


      coolneq=10**rLin

      bbccool=bbccool+coolneq*denH*dene*1e-23*Zmet

       return
       end
c========================================================================BBCOOL

c========================================================================HUNTER
      subroutine hunter(xarr,n,x,j)
c========================================================================
c  Returns index j such that X is between xarr(j) and xarr(j+1)
c  If x is outside array, j=0 or n
      dimension xarr(n)
      jbot=0
      jtop=n+1

 100  if ( (jtop-jbot).gt. 1) then 
         jmid=(jtop+jbot)/2
         if ((xarr(n).gt.xarr(1)).eqv.(x.gt.xarr(jmid))) then 
            jbot=jmid
         else 
            jtop=jmid
         endif
         goto 100
      endif
      j=jbot
      return
      end
c========================================================================HUNTER
c========================================================================Y
      function yelect(y)
c========================================================================
c      include "bbc.h"
#include "bbc.h"
      dimension y(*)

      yelect=0.
      do 10 i=1,NIONS
 10      yelect=yelect+y(i)*ncharge(i)
      return
      end
c========================================================================Y
#endif
