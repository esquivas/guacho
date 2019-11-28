#ifdef COOLINGBBC
c=======================================================================BBCRD
      subroutine bbcrd
c=======================================================================
c      INPUT:  None, except for 9 data files (see README.bbc for details)
c      OUTPUT: Common blocks of atomic data in bbc.h
c      Written by Bob Benjamin 20 Jan 2001
c------------------------------------------------------------------
#include "bbc.h"

c     -----species_rd info----------------
      character*4 label_
      character*80 explan_text

      character*6 roman(30)
      data roman/'I     ','II    ','III   ','IV    ','V     ',
     >           'VI    ','VII   ','VIII  ','IX    ','X     ',
     >           'XI    ','XII   ','XIII  ','XIV   ','XV    ',
     >           'XVI   ','XVII  ','XVIII ','XIX   ','XX    ',
     >           'XXI   ','XXII  ','XXIII ','XXIV  ','XXV   ',
     >           'XXVI  ','XXVII ','XXVIII','XXIX  ','XXX   '/


c     -----rates_rd info----------------
      character*80 labeljunk

c     -----coolneq_rd info----------------
      dimension b(NTneq), c(NTneq)
      dimension bb(NTneq-1), cc(NTneq-1)
      dimension dd(182)

c     bbc0=species  bbc1=rates(recom/ioniz) bbc2=cooling(H/He)  
c     bbc3=eq (equilibrium cooling/zbar) bbc4=D (D(Z,T))
c     bbc5=L (L(Z,T))  bbc6=Zi (ioninzing near equilibrium)
c     bbc7=Zr (recombining near equilibrium)
c     bbc8=f (ionic fractions as a function of Z,T)
      open(unit=22,file=
     & 'cool_libbbc0.table',status='unknown')
      open(unit=23,file=
     & 'cool_libbbc1.table',status='unknown')
      open(unit=24,file=
     & 'cool_libbbc2.table',status='unknown')
      open(unit=16,file=
     & 'cool_libbbc3.table',status='unknown')
      open(unit=17,file=
     & 'cool_libbbc4.table',status='unknown')
      open(unit=18,file=
     & 'cool_libbbc5.table',status='unknown')
      open(unit=19,file=
     & 'cool_libbbc6.table',status='unknown')
      open(unit=20,file=
     & 'cool_libbbc7.table',status='unknown')
      open(unit=21,file=
     & 'cool_libbbc8.table',status='unknown')
c     I changed the first three unit numbers to avoid conflict w/ out.F
c     ------------------------------------------------
c      call species_rd(10)
c     ------------------------------------------------
c     Reads in a file that contains generic information on atomic data
c     There is very little error checking. Beware.
c     INPUT :            atomic.inp
c     FORMAT OF INPUT FILE:  
c     Line 1   :   Text
c     Line 2   :   number of elements, number of molecules
c     Line 3   :   Text
c     Line 4-x :   x=number of elements+number of molecules
c                  1) Label    (4 characters max)
c                  2) Atomic # (0 for molecules)
c                  3) Charge   (0 for atom)
c                  4) Mass     (in atomic units)
c                  5) Gam      (=Cp/Cv; 5/3 for ideal monatomic gas)
c                              (        7/5 for diatomic molecules )
c                  6) Abundance(=1 for molecules                   )
c     OUTPUT:  fills up arrays in common block atdat1
c               
c-----------------------------------------------------------------------
      !print *, 'reading bbc0.table... (species info)'
      nunit=22
      
      read(nunit,1001) explan_text
1001  format(a80)
      read(nunit,*) Nel_,Nother_

      if ( (Nel_.ne.NEL).or.(Nother_.ne.NOTHER) ) then
        print *, 'Bad Nel_ or Nother_'
        stop
      endif

c     ------------------------
c     Fill up element arrays  
c     ------------------------    
      i=1
      read(nunit,1001) explan_text
      do 10 iel=1,NEL
          read(nunit,1002) label_, no_, ncharge_,amu_, gam_, abund_
1002      format(a4,2x,i3,6x,i3,3x,f9.4,f8.4,e9.2 )
 
          indexel(iel)=i
          noel(iel)=no_

          do 20 jion=1,no_+1
              no(i)=no_
              ylab(i)=label_//roman(jion)
              ncharge(i)=jion-1
              amu(i)=amu_
              gam(i)=gam_
              abund(i)=abund_
              i=i+1
 20       continue

 10   continue

c     Fill up molecule/other arrays
      do 30 iot=1,Nother_
          read(nunit, 1002) label_, no_, ncharge_, amu_, gam_, abund_
          no(i)=no_
          ylab(i)=label_//'      '
          ncharge(i)=ncharge_
          amu(i)=amu_
          gam(i)=gam_
          abund(i)=abund_
          i=i+1
 30   continue

c     Fill up arrays for electron
      ylab(i)='e-        '
      ncharge(i)=-1
      amu(i)=5.486e-4
      gam(i)=1.6667
      abund(i)=1.00

      if (i.ne.NTOT) then 
        print *, 'i does not equal NTOT!'
        stop
      endif

c     ------------------------------------------------
c      call rates_rd(11)
c     ------------------------------------------------
      !print *, 'reading bbc1.table... (ioniz/recom rates)'
      nunit=23


        read(nunit,1001) labeljunk
        read(nunit,1001) labeljunk
        read(nunit,*) ny,nT
        if ( (ny.ne.NIONS).or.(nt.ne.NTEMP)) then
           print *, 'ny/NIONS/nt/NTEMP in error!'
        endif
        do 90 iy=1,NIONS
           read(nunit,1001) labeljunk
           read(nunit,1001) labeljunk
           read(nunit,*) errtmp
           do 95 iT=1,NTEMP
              read(nunit,*) T_(iT),rcol(iy,iT)
 95        continue
 90     continue

        do 91 iy=1,NIONS
           read(nunit,1001) labeljunk
           read(nunit,1001) labeljunk
           read(nunit,*) errtmp
           do 96 iT=1,NTEMP
              read(nunit,*) T_(iT),rrec(iy,iT)
 96        continue
 91     continue

c       Spot check to make sure that temperature bins are spaced
c       uniformly in log T
        dT1=T_(2)-T_(1)
        dTf=T_(NTEMP)-T_(NTEMP-1)
        Tcomp=abs(dT1/dTf)   
        if ( (Tcomp.gt. 1.001).or.(Tcomp.lt.0.999) ) then 
          print *, 'Problem with spacing in T'
          stop
        endif

c     ------------------------------------------------
c      call cooling_rd(13)
c     ------------------------------------------------
      !print *, 'reading bbc2.table...(partial cooling tables)'
      nunit=24
        
        read(nunit,1001) labeljunk
        read(nunit,1001) labeljunk
        read(nunit,*) ny,nT        

        if( (ny.ne.NIONS).or.(nT.ne.NTEMP)) then
          print *, 'ny NIONS nT NTEMP not correct'
        endif

         do 100 iy=1,NIONS
            read(nunit,1001) labeljunk  
            read(nunit,1001) labeljunk    
            read(nunit,*) errtmp       
            do 106 iT=1,NTEMP
               read(nunit,*) T_(iT), rcool(iy,iT)
 106        continue
 100     continue


c     ------------------------------------------------
c      call coolneq_rd(16,17,18,19,20,21)
c     ------------------------------------------------
       n0=16
       n1=17
       n2=18
       n3=19
       n4=20
       n5=21

       !print *, 'reading bbc3.table... (equilibrium zbar/cooling)'
       k=0
       do i=1,NTneq+10
        read(n0,*,end=500) Tneq(i), Zeq(i),rLeq(i) 
        k=k+1
       enddo

 500   if (k.ne.NTneq) then 
        print *, 'Not NTneq entries in eq.dat'
        stop
       endif

       !print *, 'reading bbc4/5.table... (L(Z,T) and D(Z,T))'
       read(n1,*) NTneq_, NZneq_
       if ( (NTneq_.ne.NTneq).or.(NZneq_.ne.NZneq)) then
         print *, 'NTneq/NZneq not correct in unit ', n1
         print *, NTneq_, NTNeq
         print *, NZneq_, NZneq
         stop
       end if
       read(n2,*) NTneq_, NZneq_
       if ( (NTneq_.ne.NTneq).or.(NZneq_.ne.NZneq)) then
         print *, 'NTneq/NZneq not correct in unit ', n2
         print *, NTneq_, NTNeq
         print *, NZneq_, NZneq
         stop
       end if

       read(n1,*) b
       read(n2,*) b
       read(n1,*) Zneq
       read(n2,*) Zneq

       do i=1,NZneq
        read(n1,*) b
        read(n2,*) c
        do k=1,NTneq
         D(k,i)=b(k)
         rL(k,i)=c(k)
        enddo
       enddo

       !print *, 'reading bbc6/7.table... (finely gridded D(Z,T))'
       read(n3,*) NTneq_, NZneq_
       if ( (NTneq_.ne.NTneq-1).or.(NZneq_.ne.NZneqf)) then
         print *, 'NTneq/NZneqf not correct in unit ', n3
         print *, NTneq_, NTNeq
         print *, NZneq_, NZneq
         stop
       end if
       read(n4,*) NTneq_, NZneq_
       if ( (NTneq_.ne.NTneq-1).or.(NZneq_.ne.NZneqf)) then
         print *, 'NTneq/NZneqf not correct in unit ', n4
         print *, NTneq_, NTNeq
         print *, NZneq_, NZneq
         stop
       end if

       read(n3,*) bb
       read(n4,*) bb
       read(n3,*) Zf
       read(n4,*) Zf

       do i=1,NZneqf
        read(n3,*) bb
        read(n4,*) cc
        do k=1,NTneq-1
         Di(k,i)=bb(k)
         Dr(k,i)=cc(k)
        enddo
       enddo

c      -------------------------------
c      No error checking on ionization
c      tables....
c      ------------------------------
       !print *, 'reading bbc8.table... (ionization fractions)'
       read(n5,*) Nio_, NTneq_, NZneq_
       read(n5,*) dd
       read(n5,*) bb
       read(n5,*) Zneq
       do i=1,NZneq
        do k=1,NTneq
         read(n5,*) dd
         do j=1,182
          fZT(j,k,i)=dd(j)
         enddo
        enddo
       enddo


c     --------------------------------------------------
c     Close all input files...
c     --------------------------------------------------

      close(22)
      close(23)
      close(24)
      close(16)
      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      return
      end
c========================================================================BBCRD
#endif

