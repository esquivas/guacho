!     Define the number of elements, the number of other species of
!     interest, and the total number of species (where Ntot=all
!     atomic ions + molecules + electrons)
!     ----------------------------------
!     Array size parameters
!        NEL= no of elements
!        NIONS=total # of ions
!        NOTHER=# of non-ionic species
!        NTEMP=# of temperature zones
!              in atomic data tables
!     -----------------------------------
      parameter(NEL=2)
      parameter(NIONS=5)
      parameter(NOTHER=1)
      parameter(NTOT=7)
      parameter(NTEMP=250)

!     ----------------------------------
!     Atomic and mathematical constants
!     ----------------------------------
!     Grevesse and Anders abundances...     

      parameter(RHO0=  2.363848191418787E-024)  
      parameter (RGAS=8.3143e7)
      parameter (PC1=3.09e18)
      parameter (YR1=3.15e7)

!     ----------------------------
!     Tab character (ASCII(9)=^I)
!     ----------------------------
      character tab
      common /charset/tab

!     ----------------------------
!     Cooling info
!     ----------------------------
      parameter (NTneq=41) 
      parameter (NZneq=250)
      parameter (NZneqf=200)

      common /cnoneq1/Tneq(NTneq),rLeq(NTneq),Zeq(NTneq),               &
     &              Zneq(NZneq),                                        &
     &              D(NTneq,NZneq),                                     &
     &             rL(NTneq,NZneq)                                      &
      common/cnoneq2/Zf(NZneqf),Dr(NTneq-1,NZneqf),Di(NTneq-1,NZneqf)
      common/cnoneq3/fZT(182,NTneq, NZneq)
      

      character*11 ylab
!     --------------------------------
!     Atomic data arrays
!     --------------------------------
      common /atdat0/  T_(NTEMP)
      common /atdat1/  no(NTOT),ncharge(NTOT),amu(NTOT),                &
     &     gam(NTOT), abund(NTOT),Eioz(NTOT) 
      common /atdat2/  rcol(NIONS,NTEMP),rrec(NIONS,NTEMP)
      common /atdat3/  rcool(NIONS,NTEMP)
      common /atdatl/  ylab(NTOT)
      common  /eldat1/ indexel(NEL),noel(NEL)      

